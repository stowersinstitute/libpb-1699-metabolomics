from cavefinomics import AstyanaxMe
from argparse import ArgumentParser
from numpy import log10, linspace
from itertools import repeat, chain
from seaborn import catplot, FacetGrid, palplot, color_palette, set_palette
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import matplotlib.ticker as plticker
from itertools import chain
flatten = chain.from_iterable
import os
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

set_palette('Set2')

parser = ArgumentParser(description="Plot metabolites of interest.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--output", type=str, help="Output.")
parser.add_argument("--output-dir-extra", type=str, help="Output for smaller plots.")
args = parser.parse_args()

ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    sample_sheet_path=args.sample_sheet,
    hmdb_file=args.hmdb,
    )

compounds = [
    'ascorbic acid',
    'dehydroascorbic acid',
    'glutathione',
    'alpha-ketoglutarate',
    'nicotinamide',
    #'nicotinic acid',
    'orotic acid',
    #'phosphoethanolamine',
  ]

if args.exclude_outlier:
    outlier = 'without-outliers'
else:
    outlier = 'with-outliers'
pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = ['4d Starved', '30d Starved', 'Refed']
conditions_really_short = ['4', r'30', 'R']
comparisons = {'PvS':('Pachon','Surface'), 'TvS':('Tinaja','Surface'), 'PvT':('Pachon','Tinaja')}
condmap = {'30d':'30d Starved', '4d':'4d Starved', 'Ref': 'Refed'}
categories = {"Aminoacids":'Amino acids',"Carbohydrates_-CCM": 'Carbohydrates / CCM',"Fattyacids":'Fatty acids',"Misc._-_sec.metabolites":'Misc',"Nucleotides":'Nucleotides'}
datasets = []
sig = {}
up = {}
for cat in categories:
    for tissue in tissues:
        for cond in condmap:
            for comp in comparisons:
                d = read_csv(f'out/work/primary/glm/singlefactor/{outlier}/{cat}/{tissue}/{cond}/{comp}.csv',index_col=0)
                d['Category'] = cat
                d['Tissue'] = tissue
                d['Condition'] = cond
                d['Comparison'] = comp
                for m in d.index:
                    if d.loc[m,'Pr(>|z|)'] < 0.05:
                        sig[m,tissue,cond,comp] = True
                    else:
                        sig[m,tissue,cond,comp] = False
                    if d.loc[m,'Estimate'] > 0.:
                        up[m,tissue,cond,comp] = True
                    else:
                        up[m,tissue,cond,comp] = False
                datasets.append(d)
significance_data = concat(datasets,axis=0).dropna()


astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
astyanax_data.columns = (' '.join((c,str(n))) for c,n in zip(astyanax_data.columns,chain.from_iterable(repeat(range(1,6+1),9*3))))
astyanax_data = astyanax_data.rename(ame.get_kegg_to_name_map(), axis=0)
astyanax_data = astyanax_data.loc[compounds]

outliers = ['Tinaja Liver Refed 6']

def process_outlier(exclude,subset):
    for outlier in outliers:
        if exclude and outlier in subset.columns:
            subset = subset.loc[:,~subset.columns.str.contains(outlier)]
    return subset

astyanax_data = process_outlier(args.exclude_outlier,astyanax_data)

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = ['30d Starved', '4d Starved', 'Refed']
condmap = {v:k for k,v in {'30d':'30d Starved', '4d':'4d Starved', 'Ref': 'Refed'}.items()}
comparisons = ['PvS','TvS','PvT']

data = []
for pop in pops:
    for tissue in tissues:
        for condition in conditions:
            for compound in compounds:
                subset = astyanax_data.loc[compound,
                  astyanax_data.columns.str.contains(pop) & astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition)]
                for val in subset:
                    data.append({'Population':pop if pop != 'Pachon' else 'Pachón','Tissue':tissue,'Condition':condmap[condition],'Compound':compound,'Value':val})
data = DataFrame(data)

# plot all compounds at once
catplot('Condition', 'Value', data=data, kind='point', row='Compound',row_order=compounds,hue='Population',col='Tissue',sharey='row',palette=['firebrick','goldenrod','dodgerblue'],capsize=0.1,height=4.)

def fix_catplot(cpds):
    for ax,compound in zip(plt.gcf().get_axes(),flatten([c]*3 for c in cpds)):
        if data[data['Compound'] == compound]['Value'].max() > 1e6:
            ax.yaxis.set_major_locator(plticker.MultipleLocator(1e6))
        elif data[data['Compound'] == compound]['Value'].max() > 1e5:
            ax.yaxis.set_major_locator(plticker.MultipleLocator(1e5))
        elif data[data['Compound'] == compound]['Value'].max() > 1e4:
            ax.yaxis.set_major_locator(plticker.MultipleLocator(1e4))
        else:
            ax.yaxis.set_major_locator(plticker.MultipleLocator(2e3))
        ax.yaxis.set_major_formatter(EngFormatter(unit="", places=0, sep="\N{THIN SPACE}"))
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.tick_params(axis='x', labelsize=14.)
        ax.tick_params(axis='y', labelsize=16.)
        #https://cduvallet.github.io/posts/2018/11/facetgrid-ylabel-access
    for i in range(len(cpds)):
        plt.gcf().get_axes()[i*3].set_ylabel(cpds[i], fontsize=20, fontweight='bold')
    for ax,name in zip(plt.gcf().get_axes(),chain.from_iterable((tissues,['']*3*(len(cpds)-1)))):
        ax.set_title(name, fontsize=24, fontweight='bold')
    for ax in plt.gcf().get_axes():
        ax.set_xticklabels(ax.get_xticklabels(),fontsize=18,fontweight='bold')
    for i,cpd in enumerate(cpds):
        for j,tissue in enumerate(tissues):
            for x,cond in enumerate(condmap.values()):
                stars = []
                for comp,color in zip(comparisons,['firebrick','goldenrod','#35d0ba']):
                    pval = significance_data.loc[(significance_data['Tissue'] == tissue) & (significance_data['Condition'] == cond) & (significance_data['Comparison'] == comp)]['Pr(>|z|)'][cpd]
                    if pval < 0.05:
                        stars.append(color)
                w = 0.1
                for star_color,offset in zip(stars,linspace(-len(stars)*w/2.,len(stars)*w/2.,len(stars))):
                        plt.gcf().get_axes()[i*3+j].text(x+offset,0.985,'*',size=18,fontweight='bold',ha='center',va='center',color=star_color,transform=plt.gcf().get_axes()[i*3+j].get_xaxis_transform())

fix_catplot(compounds)

plt.savefig(args.output)

# make smaller plots for tiling
if args.output_dir_extra is not None:
    compounds = [
      ['ascorbic acid','dehydroascorbic acid'],
      ['glutathione','alpha-ketoglutarate'],
      ['nicotinamide','orotic acid'],
      ]
    for k,cpds in enumerate(compounds,start=1):
        d = data[data['Compound'].isin(cpds)]
        catplot('Condition', 'Value', data=d, kind='point', row='Compound',row_order=cpds,hue='Population',col='Tissue',sharey='row',palette=['firebrick','goldenrod','dodgerblue'],capsize=0.1,height=4.,legend=None)

        fix_catplot(cpds)

        handles, labels = plt.gcf().get_axes()[0].get_legend_handles_labels()

        plt.savefig(os.path.join(args.output_dir_extra,f'metabolites-of-interest{k}.pdf'),bbox_inches='tight',transparent=True,pad_inches=0)

# plot just legend
fig,ax = plt.subplots(figsize=(2.5, 2.5))
ax.legend(handles, labels, loc='center')
# hide graphics
ax.set_yticks([], minor=[])
ax.set_xticks([], minor=[])
ax.patch.set_visible(False)
for s in ["top", "bottom", "left", "right"]:
    ax.spines[s].set_visible(False)
plt.savefig(os.path.join(args.output_dir_extra,f'metabolites-of-interest-legend.pdf'),bbox_inches='tight',transparent=True,pad_inches=0)
