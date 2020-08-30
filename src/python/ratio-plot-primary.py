from argparse import ArgumentParser
from pandas import read_csv
from cavefinomics import AstyanaxMe
from seaborn import heatmap, distplot
from seaborn import catplot, FacetGrid
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import matplotlib
from numpy import log10, isfinite, any, all, quantile, nan, isnan
from pandas import DataFrame, concat, melt
import pathlib
import os
from cavefinomics import run_opls
from sklearn.preprocessing import StandardScaler, Normalizer
from numpy import log10, mean, isnan, quantile, nan, array
import matplotlib.patches as patches
import json
#plt.style.use('seaborn')
#plt.style.use('Solarize_Light2')
#plt.style.use('fivethirtyeight')
#matplotlib.rc('image', cmap='Pastel2')
#plt.set_cmap('jet')
#plt.set_cmap('Pastel2')
#https://colorhunt.co/palette/192017
#matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#99b898", "#feceab", "#ff847c","#e84a5f"])
#https://colorhunt.co/palette/191947
#matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#111d5e", "#c70039", "#f37121","#ffbd69"])
#https://colorhunt.co/palette/192164
#matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#e7305b", "#e2979c", "#f7f5dd","#9bdeac"])
#https://colorhunt.co/palette/189676
#matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#726a95", "#709fb0", "#a0c1b8","#f4ebc1"])
#https://colorhunt.co/palette/201413
#https://colorhunt.co/palette/201843
#https://colorhunt.co/palette/195532
#https://colorhunt.co/palette/179483
#https://colorhunt.co/palette/178369

parser = ArgumentParser(description="PCA vs mammals.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--output", type=str, help="Output file.")
args = parser.parse_args()


ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    sample_sheet_path=args.sample_sheet,
    hmdb_file=args.hmdb,
    )


# get data and compute common metabolites
astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
weights_series = ame.column_table['Mass (mg)']
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
weights_series.index = [' '.join(u) for u in ame.treatment_descriptors]
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
weights_series = weights_series.loc[['pools' not in c for c in weights_series.index]]
kegg_to_category = {kegg:cat for cat,keggs in ame.compounds_by_category_from_dataset.items() for kegg in keggs}
#print(astyanax_data)
astyanax_data['Category'] = astyanax_data.apply(lambda u: kegg_to_category[u.name],axis=1)
print(astyanax_data)

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
polarities = ['positive','negative']
conditions = ['4d Starved', '30d Starved', 'Refed']
conditions_short = ['4d\nStarved', '30d\nStarved', 'Refed']

outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def process_outlier(subset):
    for o in outliers:
        subset.loc[:,subset.columns.str.contains(o)] = nan
    return subset

catdata = []
#classs = []
sum_by_category = astyanax_data.groupby('Category').sum()

for tissue in tissues:
    for i,condition in zip(range(3),conditions):
        for j,pop in zip(range(1,4),pops):
            subset = sum_by_category.loc[:,sum_by_category.columns.str.contains(tissue) & sum_by_category.columns.str.contains(condition) & sum_by_category.columns.str.contains(pop)]

            def process_df(d):
                #d.index.name = 'Category'
                d = d.reset_index()
                d.columns = ['Category'] + ['Intensity']*(len(d.columns)-1)
                d.insert(0,'Tissue',[tissue]*len(d.index))
                d.insert(1,'Population',[pop]*len(d.index))
                d.insert(2,'Condition',[condition]*len(d.index))
                #d['Category'] = d['Category'].apply(lambda u: u.split('[')[0].strip())
                #print(d)
                return d

            catdata.append(process_df(subset))

catdata = concat(catdata,axis=0)
#classs = concat(classs,axis=0)
catdata = catdata.groupby(['Tissue','Category','Population','Condition']).sum().reset_index()
#print(catdata)
#print(cats)
#print(classs)

cattypes = ['Categories','Classes']


#class_subsets = ['Saturated Fatty Acids', 'Monounsaturated Fatty Acids', 'Polyunsaturated Fatty Acids']
#class_subsets = ['Ceramides','Fatty Acids and Conjugates','Glycerophosphocholines','Glycerophosphoethanolamines','Neutral glycosphingolipids','Sphingoid bases','Triradylglycerols']
#class_subsets = ['Ceramides','Fatty Acids and Conjugates','Glycerophosphocholines','Glycerophosphoethanolamines','Sphingoid bases','Triradylglycerols']
#class_subsets = ['Ceramides','Fatty Acids and Conjugates','Glycerophosphocholines','Glycerophosphoethanolamines','Triradylglycerols']
class_subsets = list(ame.compounds_by_category_from_dataset.keys())
class_renamer = {
  'Fatty Acids and Conjugates': 'Fatty Acids /\nConjugates',
  'Glycerophosphoinositols': 'Glycerophospho-\ninositols',
  'Glycerophosphoethanolamines': 'Glycerophospho-\nethanolamines',
}

data = {}
flattened = {}
for tissue in tissues:
    for condition in conditions:
        subset = {}
        for pop in pops:
            df = catdata.loc[(catdata['Tissue'] == tissue) & (catdata['Population'] == pop) & (catdata['Condition'] == condition)]
            df = df.set_index('Category')
            df = df.iloc[:,3:]
            df.columns = [pop]*len(df.columns)
            subset[pop] = df

        df = concat([subset[pop].transpose() for pop in pops],axis=0)
        df.index.name = 'Population'

        def flatten(d):
            d = d.stack().reset_index()
            cols = list(d.columns)
            cols[-1] = 'LogIntensity'
            d.columns = cols
            d.insert(1, 'Condition', condition)
            d.insert(2, 'Tissue', tissue)
            return d
        flattened[tissue,condition] = flatten(df)
        df = df.reset_index()
        df.insert(1, 'Condition', condition)
        df.insert(2, 'Tissue', tissue)

        #print(tissue,condition)
        data[tissue,condition] = df

data = concat((data[tissue,condition] for tissue in tissues for condition in conditions),axis=0)
flattened = concat((flattened[tissue,condition] for tissue in tissues for condition in conditions),axis=0)


def make_fig(class_subsets, class_renamer, name):
    #print(class_subsets)
    gridspec_kw = {"height_ratios":[1.]*3, "width_ratios" : [1.,3.,3.,3.]}
    fig,ax = plt.subplots(nrows=3,ncols=len(tissues)+1,sharex=False,sharey=False,gridspec_kw=gridspec_kw,figsize=(12.,12.))
    #fig.suptitle(f'Lipid Composition')
    c = {}
    for i,tissue in zip(range(3),tissues):
        # plot by categories
        #ax[k*4,2].text(0.5,0.5,tissue,size=12,ha='center',va='center',transform=ax[k*4,2].transAxes)
        ax[i,0].text(0.5,0.5,tissue,size=12,rotation=90.,ha='center',va='center',transform=ax[i,0].transAxes)
        # hide graphics
        ax[i,0].set_yticks([], minor=[])
        ax[i,0].set_xticks([], minor=[])
        ax[i,0].patch.set_visible(False)
        for s in ["top", "bottom", "left", "right"]:
            ax[i,0].spines[s].set_visible(False)
        for j,pop in zip(range(1,4),pops):
            d = data.loc[(data['Tissue'] == tissue) & (data['Population'] == pop)]
            d = d.groupby('Condition').mean()
            #print(d)
            #print(class_subsets)
            #condition_totals = d[class_subsets].sum(axis=1)
            condition_totals = d.sum(axis=1)
            #print(condition_totals)
            #stop
            d = d.rename(class_renamer,axis=1)
            #print(d)
            #print(condition_totals)
            for cls in class_subsets:
                if cls in class_renamer:
                    cls = class_renamer[cls]
                #print(cls)
                c[tissue,pop,cls] = d[cls]/condition_totals
                #d = d[class_subsets].mean()
                #print(tissue,pop,cls)
                if not all(isfinite(c[tissue,pop,cls].values)):
                    c[tissue,pop,cls] = [0.]*len(c)
                    #print(tissue,condition,pop,'Non finite values')
                    #continue
                #https://python-graph-gallery.com/12-stacked-barplot-with-matplotlib/
                #print(c)
    for i,tissue in zip(range(3),tissues):
        for j,condition in zip(range(1,4),conditions):
            last_bar = None
            for cls in class_subsets:
                if cls in class_renamer:
                    cls = class_renamer[cls]
                try:
                    vals = {pop: c[tissue,pop,cls][condition] for pop in pops}
                except:
                    continue
                kwds = {}
                if i==0 and j==1:
                    kwds['label'] = cls
                if last_bar is not None:
                    kwds['bottom']=array(last_bar)
                    last_bar += array(list(vals.values()))
                else:
                    last_bar = array(list(vals.values()))
                #for _,cv in zip(range(3),c):
                #print(list(vals.keys()))
                #print(list(vals.values()))
                ax[i,j].bar(list(vals.keys()), list(vals.values()), tick_label=pops, **kwds) #
                if i == 0:
                    ax[i,j].set_title(conditions_short[j-1])
            #ax[i,j].bar(pops, array([1.]*3)-last_bar, tick_label=pops, bottom=last_bar, label='Other' if i==0 and j==1 else "") #
            # hide graphics
            ax[i,j].set_yticks([], minor=[])
            #ax[i,j].set_xticks([], minor=[])
            ax[i,j].patch.set_visible(False)
            for s in ["top", "bottom", "left", "right"]:
                ax[i,j].spines[s].set_visible(False)

    fig.legend()

    #pathlib.Path('/tmp/classes/bar').mkdir(parents=True, exist_ok=True)
    plt.savefig(args.output,bbox_inches='tight',transparent=True,pad_inches=0)

#make_fig(list(data.columns)[3:], {}, 'all')
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#363062", "#4d4c7d", "#827397","#d8b9c3","#ffb6b6","#fde2e2","#aacfcf","#679b9b","#318fb5"])
make_fig(class_subsets, class_renamer, 'subset')
