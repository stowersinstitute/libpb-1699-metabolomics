from cavefinomics import AstyanaxMe
from argparse import ArgumentParser
from numpy import log10
from itertools import repeat, chain
from seaborn import (catplot, FacetGrid, palplot,
                     color_palette, set_palette, heatmap)
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import matplotlib.ticker as plticker
from pandas import concat

SMALL_SIZE = 12
MEDIUM_SIZE = 12
BIGGER_SIZE = 12

#https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
#colors = ['firebrick','goldenrod','mediumseagreen','rosybrown']
#https://colorhunt.co/palette/184190
#colors = ['#0e9aa7','#3da4ab','#f6cd61','#fe8a71']
#https://colorhunt.co/palette/189889
colors = ['#d92027','#ff9234','#ffcd3c','#35d0ba']
#https://colorhunt.co/palette/183823
#colors = ['#ffc4a3','#ff9a76','#f96d80','#bb596b']

#palplot(color_palette('Set1'))
set_palette('Set2')

parser = ArgumentParser(description="PCA vs mammals.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--level", type=float, help="P value cutoff.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--output", type=str, help="Output.")
args = parser.parse_args()

ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    sample_sheet_path=args.sample_sheet,
    hmdb_file=args.hmdb,
    )

compounds = [
    'glucose-1-phosphate',
    'glucose-6-phosphate',
    'galactose-6-phosphate',
    'fructose-1-phosphate',
    'fructose-6-phosphate',
    'ribose-5-phosphate',
    'ribulose-5-phosphate',
    'glucose',
    'fructose',
    'glucuronic acid',
    'gluconic acid',
  ]

astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
astyanax_data.columns = (' '.join((c,str(n))) for c,n in zip(astyanax_data.columns,chain.from_iterable(repeat(range(1,6+1),9*3))))
astyanax_data = astyanax_data.rename(ame.get_kegg_to_name_map(), axis=0)
astyanax_data = astyanax_data.loc[compounds]

outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def process_outlier(exclude,subset):
    for outlier in outliers:
        if exclude and outlier in subset.columns:
            subset = subset.loc[:,~subset.columns.str.contains(outlier)]
    return subset

astyanax_data = process_outlier(args.exclude_outlier,astyanax_data)

pops = ['Pachon', 'Tinaja', 'Surface']
pops2 = ['Pachón', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
#conditions = ['30d Starved', '4d Starved', 'Refed']
conditions = {'4d':'4d Starved', '30d':'30d Starved', 'Ref':'Refed'}
condmap = {v:k for k,v in {'30d':'30d Starved', '4d':'4d Starved', 'Ref': 'Refed'}.items()}
comparisons = ['PvS','TvS','PvT']
categories = {"Aminoacids":'Amino acids',"Carbohydrates_-CCM": 'Carbohydrates / CCM',"Fattyacids":'Fatty acids',"Misc._-_sec.metabolites":'Misc',"Nucleotides":'Nucleotides'}

data = []
for pop in pops:
    for tissue in tissues:
        for condition in conditions.values():
            for compound in compounds:
                subset = astyanax_data.loc[compound,
                  astyanax_data.columns.str.contains(pop) & astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition)]
                for val in subset:
                    data.append({'Population':pop if pop != 'Pachon' else 'Pachón','Tissue':tissue,'Condition':condmap[condition],'Compound':compound,'Value':val})
data = DataFrame(data)

# significance data
datasets = []
sig = {}
up = {}
if args.exclude_outlier:
    outlier_text = 'no-outliers'
else:
    outlier_text = 'outliers'
for cat in categories:
    for tissue in tissues:
        for cond in conditions:
            for comp in comparisons:
                d = read_csv(f'out/work/primary/glm/singlefactor/{outlier_text}/{cat}/{tissue}/{cond}/{comp}.csv',index_col=0)
                d['Category'] = cat
                d['Tissue'] = tissue
                d['Condition'] = cond
                d['Comparison'] = comp
                for m in d.index:
                    if d.loc[m,'Pr(>|z|)'] < args.level:
                        sig[m,tissue,cond,comp] = True
                    else:
                        sig[m,tissue,cond,comp] = False
                    if d.loc[m,'Estimate'] > 0.:
                        up[m,tissue,cond,comp] = True
                    else:
                        up[m,tissue,cond,comp] = False
                datasets.append(d)
sig_data = concat(datasets,axis=0).dropna()
print(sig_data)


fig,ax = plt.subplots(nrows=len(set(data['Compound'])),ncols=len(tissues)*len(conditions),figsize=(12.,12.))

#print(data)
for tissue in tissues:
    for cond_label,condition in conditions.items():
        data2d = data.loc[(data['Tissue'] == tissue) & (data['Condition'] == cond_label)]
        data2d = data2d.drop('Tissue',axis=1).drop('Condition',axis=1)
        data2d = data2d.pivot_table(index='Compound',columns='Population')
        #https://stackoverflow.com/questions/35678874/normalize-rows-of-pandas-data-frame-by-their-sums/35679163
        data2d = data2d.div(data2d.max(axis=1),axis=0)
        #https://stackoverflow.com/questions/39273441/flatten-pandas-pivot-table
        data2d.columns = data2d.columns.to_series().str.join('_').str.replace('Value_','')
        #print(list(data2d.columns))
        data2d = data2d[pops2]
        data2d = data2d.loc[compounds,:]
        print(data2d)
        print(len(data2d.index))

        #print(list(sig_data.columns))
        #print(sig_data['Comparison'] == 'PvS')
        sig2d = sig_data.loc[(sig_data['Comparison'] == 'PvS') & (sig_data['Condition'] == cond_label) & (sig_data['Tissue'] == tissue)].loc[data2d.index,'Pr(>|z|)']
        sig2d = sig2d.loc[compounds]
        print(sig2d)
        #print(data2d)
        stop
        #heatmap(data2d, vmin=0., vmax=1., annot=array([annot]), fmt = '', ax=axs[j,i], cbar=False, xticklabels=j==len(categories)-1, yticklabels=i==0, cmap='coolwarm')



plt.savefig(args.output)
