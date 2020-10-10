from argparse import ArgumentParser
import os

from numpy import array, log10
from pandas import DataFrame, Series
from sklearn.preprocessing import StandardScaler, Normalizer
import matplotlib.pyplot as plt
import matplotlib
from seaborn import heatmap
from sklearn.decomposition import PCA as PCA
from functools import reduce
from scipy.stats import chisquare, fisher_exact
import math
from pprint import pprint
from collections import defaultdict

#from cavefinomics import AstyanaxMe, MaMammalianMe
from cavefinomics import AstyanaxLi

parser = ArgumentParser(description="PCA for lipids.")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument("--lipidmaps-json", type=str, help="Lipidmaps JSON.")
parser.add_argument("--exclude-outlier", type=int, help="Exclude the single Tinaja outlier?")
parser.add_argument("--output", type=str, help="Output lipids PCA.")
args = parser.parse_args()

if not args.exclude_outlier:
    raise RuntimeError('should exclude outlier')

short_groups = ['4d Starved', 'Refed']
long_groups = ['30d Starved', 'Refed']
exclude_outlier = args.exclude_outlier

ali = AstyanaxLi(
    lipids_normalized=args.lipids_normalized,
    lipidmaps_js=args.lipidmaps_json,
    )

# get data and compute common metabolites
astyanax_data = ali.lmdata.iloc[:,ali.non_numeric_cols-1:]
astyanax_data = astyanax_data.apply(log10)


pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = {'30d':'30d Starved', '4d':'4d Starved', 'Ref':'Refed'}
comparisons = {'PvS':('Pachon','Surface'),'TvS':('Tinaja','Surface'),'PvT':('Pachon','Tinaja')}
cattypes = {'Class':'class', 'Category':'category'}

cattype = 'Category'

categories = list(sorted(set(ali.lmdata[cattype])))
outliers = ['Tinaja Liver Refed 6']

colors = {
  'Pachon': 'firebrick',
  'Tinaja': 'goldenrod',
  'Surface': 'dodgerblue'
}

#https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

# already added in lipid data
def add_replicate_nums(subset):
    subset.columns = (f'{c} {n}' for n,c in zip(list(range(1,7))*9,subset.columns))
    return subset

def fix_cols(df,tissue):
    d = df.copy()
    def j(c):
        return ' '.join(c)
    #print(d.columns)
    d.columns = (f'{c[0]} {tissue} {j(c[1:-2])} {n}' for c,n in zip((cc.split() for cc in d.columns),list(range(1,7))*27))
    d = d.reindex(sorted(d.columns), axis=1)
    #print(d.columns)
    return d

if True:
    height_ratios = [1.]*len(categories)
else:
    height_ratios = [float(len(ame.compounds_by_category_from_dataset[c])) for c in categories]
gridspec_kw = {"height_ratios":height_ratios, "width_ratios" : [3.,3.,3.,1.,1]}

def process_outlier(exclude,subset):
    for outlier in outliers:
        if exclude and outlier in subset.columns:
            subset = subset.loc[:,~subset.columns.str.contains(outlier)]
    return subset

remap_cond = {
  '4d Starved': '4d Fasted',
  '30d Starved': '30d Fasted',
  }

fig,ax = plt.subplots(nrows=len(categories),ncols=5,figsize=(12, 8), gridspec_kw=gridspec_kw)
fig.suptitle('Lipids',fontsize=24,fontweight='bold',y=1.0)
for i,category in zip(range(len(categories)),categories):
    for j,tissue in zip(range(3),['Brain', 'Muscle', 'Liver']):
        # get subset
        subset = process_outlier(args.exclude_outlier,fix_cols(DataFrame(astyanax_data.loc[
            ali.lmdata[cattype] == category,
            astyanax_data.columns.str.contains(tissue)]),tissue))

        pca_data = StandardScaler().fit_transform(subset.transpose())
        pca = PCA(n_components=2)
        pca_analysis = pca.fit(pca_data)
        tf_data = DataFrame(pca_analysis.transform(pca_data), columns=['C1','C2'], index=subset.columns)

        for pop in ['Surface', 'Tinaja', 'Pachon']:
            for color,feeding_state in zip([adjust_lightness(colors[pop],c) for c in [0.5, 1.0, 1.5]], ['4d Starved', '30d Starved', 'Refed']):
                p = tf_data.iloc[subset.columns.str.contains(feeding_state) & subset.columns.str.contains(pop)]
                #print(pop,feeding_state,tissue,p)
                label = ', '.join((pop, remap_cond[feeding_state] if feeding_state in remap_cond else feeding_state))
                ax[i,j].scatter(p['C1'],p['C2'], label=label, color=color)
        handles, labels = ax[i,j].get_legend_handles_labels()

        if i == 0:
            #https://stackoverflow.com/questions/12444716/how-do-i-set-the-figure-title-and-axes-labels-font-size-in-matplotlib
            ax[i,j].set_title(tissue,fontsize=24,y=1.15)

for i,category in zip(range(len(categories)),categories):
    for j in range(3):
        ax[i,j].set_yticks([], minor=[])
        ax[i,j].set_xticks([], minor=[])
    xpos = -0.1
    if category == 'Glycerophospholipids':
        category = 'Glycerophospho-\nlipids'
    ax[i,3].text(xpos,0.5,category,size=16,rotation=0.,fontweight='bold',ha='left',va='center',transform=ax[i,3].transAxes)
    # hide graphics
    ax[i,3].set_yticks([], minor=[])
    ax[i,3].set_xticks([], minor=[])
    ax[i,3].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        ax[i,3].spines[s].set_visible(False)
    # hide graphics
    ax[i,4].set_yticks([], minor=[])
    ax[i,4].set_xticks([], minor=[])
    ax[i,4].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        ax[i,4].spines[s].set_visible(False)

#fig.legend(handles, labels, loc='center right')

plt.savefig(args.output)

