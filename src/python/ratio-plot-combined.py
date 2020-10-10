from argparse import ArgumentParser
from pandas import read_csv
from cavefinomics import AstyanaxLi, AstyanaxMe
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
from numpy import (log10, mean, isnan, quantile, nan, array, arange)
import matplotlib.patches as patches
import json
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

parser = ArgumentParser(description="Ratio of metabolite categories.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument("--lipidmaps-json", type=str, help="Lipidmaps JSON.")
parser.add_argument("--lipidmaps-fa", type=str, help="Lipidmaps fa.")
parser.add_argument("--output", type=str, help="Output file.")
parser.add_argument("--output-legend", type=str, help="Output legend file.")
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
astyanax_data['Category'] = astyanax_data.apply(lambda u: kegg_to_category[u.name],axis=1)

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
polarities = ['positive','negative']
conditions = ['4d Starved', '30d Starved', 'Refed']
conditions_short = ['4d\nStarved', '30d\nStarved', 'Refed']

outliers = ['Tinaja Liver Refed 6']

def process_outlier(subset):
    for o in outliers:
        subset.loc[:,subset.columns.str.contains(o)] = nan
    return subset

catdata = []
sum_by_category = astyanax_data.groupby('Category').sum()

for tissue in tissues:
    for i,condition in zip(range(3),conditions):
        for j,pop in zip(range(1,4),pops):
            subset = sum_by_category.loc[:,sum_by_category.columns.str.contains(tissue) & sum_by_category.columns.str.contains(condition) & sum_by_category.columns.str.contains(pop)]

            def process_df(d):
                d = d.reset_index()
                d.columns = ['Category'] + ['Intensity']*(len(d.columns)-1)
                d.insert(0,'Tissue',[tissue]*len(d.index))
                d.insert(1,'Population',[pop]*len(d.index))
                d.insert(2,'Condition',[condition]*len(d.index))
                return d

            catdata.append(process_df(subset))

catdata = concat(catdata,axis=0)
catdata = catdata.groupby(['Tissue','Category','Population','Condition']).sum().reset_index()

with open(args.lipidmaps_fa) as f:
  lmfa = json.load(f)

ali = AstyanaxLi(
  lipids_normalized=args.lipids_normalized,
  lipidmaps_js=args.lipidmaps_json,
  )

def get_lipidmap_id(inchikey):
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
    d.columns = (f'{c[0]} {tissue} {j(c[1:])} {n}' for c,n in zip((cc.split() for cc in d.columns),list(range(1,7))*27))
    if args.exclude_outlier:
        d = process_outlier(d)
    d = d.reindex(sorted(d.columns), axis=1)
    return d

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
                if mc != "Fatty Acids and Conjugates [FA01]":
                    return mc
                else:
                    fa = lmfa[lmid]
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
        d = ali.normalized[tissue,polarity]
        d = fix_cols(d,tissue)
        d['CATEGORY'] = lmd['CATEGORY']
        sum_by_category = d.groupby('CATEGORY').sum()

        d = ali.normalized[tissue,polarity]
        d['MAIN_CLASS'] = lmd['MAIN_CLASS']
        sum_by_class = d.groupby('MAIN_CLASS').sum()

        for i,condition in zip(range(3),conditions):
            for j,pop in zip(range(1,4),pops):
                category_subset = sum_by_category.loc[:,sum_by_category.columns.str.contains(condition) & sum_by_category.columns.str.contains(pop)]
                classes_subset = sum_by_class.loc[:,sum_by_class.columns.str.contains(condition) & sum_by_class.columns.str.contains(pop)]

                def process_df(d):
                    d = d.reset_index()
                    d.columns = ['Category'] + ['Intensity']*(len(d.columns)-1)
                    d.insert(0,'Tissue',[tissue]*len(d.index))
                    d.insert(1,'Polarity',[polarity]*len(d.index))
                    d.insert(2,'Population',[pop]*len(d.index))
                    d.insert(3,'Condition',[condition]*len(d.index))
                    d['Category'] = d['Category'].apply(lambda u: u.split('[')[0].strip())
                    return d

                category_subset,classes_subset = (process_df(d) for d in [category_subset,classes_subset])

                cats.append(category_subset)
                classs.append(classes_subset)

cats = concat(cats,axis=0)
classs = concat(classs,axis=0)
cats = cats.groupby(['Tissue','Category','Population','Condition']).sum().reset_index()
classs = classs.groupby(['Tissue','Category','Population','Condition']).sum().reset_index()

cattypes = ['Categories','Classes']

cattypes = ['Categories','Classes']


lipid_subsets = ['Ceramides','Fatty Acids and Conjugates','Glycerophosphocholines','Glycerophosphoethanolamines','Triradylglycerols']
primary_subsets = list(ame.compounds_by_category_from_dataset.keys())
primary_renamer = {'Fatty\nacids': 'Primary fatty acids'}
lipid_renamer = {
  'Fatty Acids and Conjugates': 'Fatty Acids /\nConjugates',
  'Glycerophosphoinositols': 'Glycerophospho-\ninositols',
  'Glycerophosphoethanolamines': 'Glycerophospho-\nethanolamines',
}

for d,c in zip([cats,classs],cattypes):
    lipid_data = {}
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

            lipid_data[tissue,condition] = df

    lipid_data = concat((lipid_data[tissue,condition] for tissue in tissues for condition in conditions),axis=0)
    lipid_flattened = concat((flattened[tissue,condition] for tissue in tissues for condition in conditions),axis=0)

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

        data[tissue,condition] = df

data = concat((data[tissue,condition] for tissue in tissues for condition in conditions),axis=0)
flattened = concat((flattened[tissue,condition] for tissue in tissues for condition in conditions),axis=0)

bar_width = 0.35
width = 0.35

remap_cond = {
  '4d\nStarved': '4d Fasted',
  '30d\nStarved': '30d Fasted',
  }

def make_fig(primary_subsets, lipid_subsets, primary_renamer, lipid_renamer, name):
    gridspec_kw = {"height_ratios":[1.]*3, "width_ratios" : [0.5,3.,3.,3.]}
    fig,ax = plt.subplots(nrows=3,ncols=len(tissues)+1,sharex=False,sharey=False,gridspec_kw=gridspec_kw,figsize=(12.,12.))
    c = {}
    lc = {}
    for i,tissue in zip(range(3),tissues):
        # plot by categories
        ax[i,0].text(0.95,0.5,tissue,size=18,fontweight='bold',rotation=90.,ha='center',va='center',transform=ax[i,0].transAxes)
        # hide graphics
        ax[i,0].set_yticks([], minor=[])
        ax[i,0].set_xticks([], minor=[])
        ax[i,0].patch.set_visible(False)
        for s in ["top", "bottom", "left", "right"]:
            ax[i,0].spines[s].set_visible(False)
        for j,pop in zip(range(1,4),pops):
            d = data.loc[(data['Tissue'] == tissue) & (data['Population'] == pop)]
            d = d.groupby('Condition').mean()
            condition_totals = d.sum(axis=1)
            d = d.rename(primary_renamer,axis=1)

            # lipids
            l = lipid_data.loc[(data['Tissue'] == tissue) & (data['Population'] == pop)]
            l = l.groupby('Condition').mean()
            lipid_condition_totals = l.sum(axis=1)
            l = l.rename(lipid_renamer,axis=1)
            for cls in primary_subsets:
                if cls in primary_renamer:
                    cls = primary_renamer[cls]
                c[tissue,pop,cls] = d[cls]/condition_totals
                if not all(isfinite(c[tissue,pop,cls].values)):
                    c[tissue,pop,cls] = [0.]*len(c)
                #https://python-graph-gallery.com/12-stacked-barplot-with-matplotlib/
            for cls in lipid_subsets:
                if cls in lipid_renamer:
                    cls = lipid_renamer[cls]
                lc[tissue,pop,cls] = l[cls]/lipid_condition_totals
                if not all(isfinite(lc[tissue,pop,cls].values)):
                    lc[tissue,pop,cls] = [0.]*len(lc)
    for i,tissue in zip(range(3),tissues):
        for j,condition in zip(range(1,4),conditions):
            primary_last_bar = None
            lipids_last_bar = None
            for cls in primary_subsets:
                if cls in primary_renamer:
                    cls = primary_renamer[cls]
                try:
                    vals = {pop: c[tissue,pop,cls][condition] for pop in pops}
                except:
                    continue
                kwds = {}
                if i==0 and j==1:
                    kwds['label'] = cls
                if primary_last_bar is not None:
                    kwds['bottom']=array(primary_last_bar)
                    primary_last_bar += array(list(vals.values()))
                else:
                    primary_last_bar = array(list(vals.values()))
                #https://matplotlib.org/3.3.1/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py
                ax[i,j].bar(array(range(3)) - width/2., list(vals.values()), bar_width, tick_label=pops, **kwds) #
                if i == 0:
                    ax[i,j].set_title(remap_cond[conditions_short[j-1]] if conditions_short[j-1] in remap_cond else conditions_short[j-1])
            for cls in lipid_subsets:
                if cls in lipid_renamer:
                    cls = lipid_renamer[cls]
                try:
                    vals = {pop: lc[tissue,pop,cls][condition] for pop in pops}
                except:
                    continue
                kwds = {}
                if i==0 and j==1:
                    kwds['label'] = cls
                if lipids_last_bar is not None:
                    kwds['bottom']=array(lipids_last_bar)
                    lipids_last_bar += array(list(vals.values()))
                else:
                    lipids_last_bar = array(list(vals.values()))
                ax[i,j].bar(array(range(3)) + width/2., list(vals.values()), bar_width, tick_label=pops, **kwds) #
                if i == 0:
                    ax[i,j].set_title(remap_cond[conditions_short[j-1]] if conditions_short[j-1] in remap_cond else conditions_short[j-1],size=16,fontweight='bold')
            ax[i,j].bar(array(range(3)) + width/2., array([1.]*3)-lipids_last_bar, bar_width, tick_label=pops, bottom=lipids_last_bar, label='Other lipids' if i==0 and j==1 else "") #
            ax[i,j].set_xticks(arange(3))
            ax[i,j].set_xticklabels(ax[i,j].get_xticklabels(),fontsize='large')
            # hide graphics
            ax[i,j].set_yticks([], minor=[])
            ax[i,j].patch.set_visible(False)
            for s in ["top", "bottom", "left", "right"]:
                ax[i,j].spines[s].set_visible(False)

    #fig.legend()
    handles, labels = ax[0,1].get_legend_handles_labels()

    plt.savefig(args.output,bbox_inches='tight',transparent=True,pad_inches=0)

    # plot just legend
    fig,ax = plt.subplots(nrows=2,figsize=(2.5, 5.))

    # hide graphics
    for a in ax:
        a.set_yticks([], minor=[])
        a.set_xticks([], minor=[])
        a.patch.set_visible(False)
        for s in ["top", "bottom", "left", "right"]:
            a.spines[s].set_visible(False)

    ax[0].legend(handles[:5], labels[:5], loc='center')
    ax[0].text(0.5,1.1,'Primary\nmetabolites',size=14,fontweight='bold',ha='center',va='center',transform=ax[0].transAxes)
    ax[1].legend(handles[5:], labels[5:], loc='center')
    ax[1].text(0.5,1.05,'Lipids',size=14,fontweight='bold',ha='center',va='center',transform=ax[1].transAxes)

    plt.savefig(args.output_legend,bbox_inches='tight',transparent=True,pad_inches=0)

matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=[
  "#28215c", "#695c88", "#98799f","#d8b9c3","#ffb6b6",
  #"#726a95",
  "#709fb0", "#a0c1b8",
  #"#f4ebc1",
  "#005086","#318fb5","#111d5e","#99b898"])
make_fig(primary_subsets, lipid_subsets, primary_renamer, lipid_renamer, 'subset')
