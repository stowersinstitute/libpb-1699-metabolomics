from cavefinomics import AstyanaxMe
from argparse import ArgumentParser
from numpy import log10
from itertools import repeat, chain
from seaborn import catplot, FacetGrid, palplot, color_palette, set_palette
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
    #'maltose',
    #'sucrose',
    #'glucose',
    #'ribulose-5-phosphate',
    #'pyruvic acid',
    #'phosphoenolpyruvate',
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
    #'orotic acid',
    #'2-hydroxyglutaric acid',
  ]

astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
#astyanax_data = astyanax_data.apply(log10)
astyanax_data.columns = (' '.join((c,str(n))) for c,n in zip(astyanax_data.columns,chain.from_iterable(repeat(range(1,6+1),9*3))))
astyanax_data = astyanax_data.rename(ame.get_kegg_to_name_map(), axis=0)
#print(list(astyanax_data.index))
astyanax_data = astyanax_data.loc[compounds]

outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def process_outlier(exclude,subset):
    for outlier in outliers:
        if exclude and outlier in subset.columns:
            subset = subset.loc[:,~subset.columns.str.contains(outlier)]
    #print(subset.columns)
    return subset

astyanax_data = process_outlier(args.exclude_outlier,astyanax_data)
#print(list(astyanax_data.columns))

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
                #print(subset)
                for val in subset:
                    data.append({'Population':pop if pop != 'Pachon' else 'Pachón','Tissue':tissue,'Condition':condmap[condition],'Compound':compound,'Value':val})
data = DataFrame(data)
#print(data)

#g = FacetGrid(data,col='Compound')
#g = (g.map(catplot, row='Population',hue='Condition',col='Tissue', y='Value'))
catplot('Condition', 'Value', data=data, kind='point', row='Compound',row_order=compounds,hue='Population',col='Tissue',sharey='row',palette=['firebrick','goldenrod','dodgerblue'],capsize=0.1,height=4.)
#catplot('Condition', 'Value', data=data, kind='point', row='Compound',row_order=compounds,hue='Population',col='Tissue',sharey=False,palette=['firebrick','goldenrod','dodgerblue'],capsize=0.1,height=4.)

for ax in plt.gcf().get_axes():
    ax.yaxis.set_major_locator(plticker.MultipleLocator(2e5))
    ax.yaxis.set_major_formatter(EngFormatter(unit="", places=0, sep="\N{THIN SPACE}"))
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='x', labelsize=20.)
    ax.tick_params(axis='y', labelsize=16.)
    ax.set_yticks([], minor=[])
    #https://cduvallet.github.io/posts/2018/11/facetgrid-ylabel-access
    #print(len(ax.texts))
    #ax.texts[0].remove()
#print(plt.gcf().get_axes().shape)
for i in range(len(compounds)):
    plt.gcf().get_axes()[i*3].set_ylabel(compounds[i], fontsize='xx-large', fontweight='bold')
#plt.gcf().get_axes()[3].set_ylabel(compounds[1], fontsize='xx-large')
for ax,name in zip(plt.gcf().get_axes(),chain.from_iterable((tissues,['']*3*(len(compounds)-1)))):
    ax.set_title(name, fontsize='xx-large', fontweight='bold')


#stop
plt.savefig(args.output)
