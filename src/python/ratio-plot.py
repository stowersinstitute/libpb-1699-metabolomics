from argparse import ArgumentParser
from pandas import read_csv
from cavefinomics import AstyanaxLi
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

parser = ArgumentParser(description="PCA vs mammals.")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument("--lipidmaps-json", type=str, help="Lipidmaps JSON.")
parser.add_argument("--lipidmaps-fa", type=str, help="Lipidmaps fa.")
parser.add_argument("--output", type=str, help="Output file.")
args = parser.parse_args()

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
polarities = ['positive','negative']
conditions = ['4d Starved', '30d Starved', 'Refed']
conditions_short = ['4d\nStarved', '30d\nStarved', 'Refed']

with open(args.lipidmaps_fa) as f:
  lmfa = json.load(f)

ali = AstyanaxLi(
  lipids_normalized=args.lipids_normalized,
  lipidmaps_js=args.lipidmaps_json,
  )

def get_lipidmap_id(inchikey):
    #return inchikey
    if inchikey in ali.lipidmaps_inchikey:
        return ali.lipidmaps_inchikey[inchikey]
    elif '-'.join(inchikey.split('-')[:2]) in ali.lipidmaps_inchikey2:
        return ali.lipidmaps_inchikey2['-'.join(inchikey.split('-')[:2])]
    elif inchikey.split('-')[0] in ali.lipidmaps_inchikey1:
        return ali.lipidmaps_inchikey1[inchikey.split('-')[0]]
    return nan

def fix_cols(df,tissue):
    d = df.copy()
    def j(c):
        return ' '.join(c)
    #print(d.columns)
    d.columns = (f'{c[0]} {tissue} {j(c[1:])} {n}' for c,n in zip((cc.split() for cc in d.columns),list(range(1,7))*27))
    d = process_outlier(d)
    d = d.reindex(sorted(d.columns), axis=1)
    #print(d.columns)
    return d

outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def process_outlier(subset):
    for o in outliers:
        subset.loc[:,subset.columns.str.contains(o)] = nan
    return subset

cats = []
classs = []
for tissue in tissues:
    for polarity in polarities:
        lmd = DataFrame(ali.normalized[tissue,polarity].apply(lambda u: get_lipidmap_id(u.name),axis=1), columns=['LMID'])
        lmd['CATEGORY'] = lmd.apply(lambda u: ali.lipidmaps[u[0]]['CATEGORY'] if isinstance(u[0],str) else nan,axis=1)
        def assignfa(u):
            if isinstance(u[0],str):
                lmid = u[0]
                mc = ali.lipidmaps[lmid]['MAIN_CLASS']
                #print(mc)
                if mc != "Fatty Acids and Conjugates [FA01]":
                    return mc
                else:
                    fa = lmfa[lmid]
                    #print(fa)
                    if 'unsat' in fa:
                        if fa['unsat'] > 1:
                            return 'Polyunsaturated Fatty Acids'
                        elif fa['unsat'] == 1:
                            return 'Monounsaturated Fatty Acids'
                    elif 'sat' in fa:
                        return 'Saturated Fatty Acids'
                    return mc
            else:
                return nan
        lmd['MAIN_CLASS'] = lmd.apply(lambda u: ali.lipidmaps[u[0]]['MAIN_CLASS'] if isinstance(u[0],str) else nan,axis=1)
        #lmd['MAIN_CLASS'] = lmd.apply(assignfa,axis=1)
        d = ali.normalized[tissue,polarity]
        d = fix_cols(d,tissue)
        d['CATEGORY'] = lmd['CATEGORY']
        sum_by_category = d.groupby('CATEGORY').sum()
        print(sum_by_category)
        stop

        d = ali.normalized[tissue,polarity]
        d['MAIN_CLASS'] = lmd['MAIN_CLASS']
        sum_by_class = d.groupby('MAIN_CLASS').sum()

        for i,condition in zip(range(3),conditions):
            for j,pop in zip(range(1,4),pops):
                category_subset = sum_by_category.loc[:,sum_by_category.columns.str.contains(condition) & sum_by_category.columns.str.contains(pop)]
                #print(sum_by_category)
                #pathlib.Path('/tmp/classes/pie').mkdir(parents=True, exist_ok=True)
                classes_subset = sum_by_class.loc[:,sum_by_class.columns.str.contains(condition) & sum_by_class.columns.str.contains(pop)]

                def process_df(d):
                    #d.index.name = 'Category'
                    d = d.reset_index()
                    d.columns = ['Category'] + ['Intensity']*(len(d.columns)-1)
                    d.insert(0,'Tissue',[tissue]*len(d.index))
                    d.insert(1,'Polarity',[polarity]*len(d.index))
                    d.insert(2,'Population',[pop]*len(d.index))
                    d.insert(3,'Condition',[condition]*len(d.index))
                    d['Category'] = d['Category'].apply(lambda u: u.split('[')[0].strip())
                    #print(d)
                    return d

                category_subset,classes_subset = (process_df(d) for d in [category_subset,classes_subset])

                cats.append(category_subset)
                classs.append(classes_subset)

cats = concat(cats,axis=0)
classs = concat(classs,axis=0)
cats = cats.groupby(['Tissue','Category','Population','Condition']).sum().reset_index()
classs = classs.groupby(['Tissue','Category','Population','Condition']).sum().reset_index()
#print(cats)
#print(classs)

cattypes = ['Categories','Classes']


#class_subsets = ['Saturated Fatty Acids', 'Monounsaturated Fatty Acids', 'Polyunsaturated Fatty Acids']
#class_subsets = ['Ceramides','Fatty Acids and Conjugates','Glycerophosphocholines','Glycerophosphoethanolamines','Neutral glycosphingolipids','Sphingoid bases','Triradylglycerols']
#class_subsets = ['Ceramides','Fatty Acids and Conjugates','Glycerophosphocholines','Glycerophosphoethanolamines','Sphingoid bases','Triradylglycerols']
class_subsets = ['Ceramides','Fatty Acids and Conjugates','Glycerophosphocholines','Glycerophosphoethanolamines','Triradylglycerols']
class_renamer = {
  'Fatty Acids and Conjugates': 'Fatty Acids /\nConjugates',
  'Glycerophosphoinositols': 'Glycerophospho-\ninositols',
  'Glycerophosphoethanolamines': 'Glycerophospho-\nethanolamines',
}

for d,c in zip([cats,classs],cattypes):
    data = {}
    flattened = {}
    for tissue in tissues:
        for condition in conditions:
            subset = {}
            for pop in pops:
                df = d.loc[(d['Tissue'] == tissue) & (d['Population'] == pop) & (d['Condition'] == condition)]
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
            ax[i,j].bar(pops, array([1.]*3)-last_bar, tick_label=pops, bottom=last_bar, label='Other' if i==0 and j==1 else "") #
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
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#726a95", "#709fb0", "#a0c1b8","#f4ebc1","#005086","#318fb5","#f7d6bf","#b0cac7"])
make_fig(class_subsets, class_renamer, 'subset')
