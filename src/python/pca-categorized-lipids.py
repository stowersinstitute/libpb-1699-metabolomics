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
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the single Tinaja outlier?")
parser.add_argument("--output", type=str, help="Output lipids PCA.")
args = parser.parse_args()

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
outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

cattype = 'Category'

categories = list(sorted(set(ali.lmdata[cattype])))

if True:
    height_ratios = [1.]*len(categories)
else:
    height_ratios = [float(len(ame.compounds_by_category_from_dataset[c])) for c in categories]
gridspec_kw = {"height_ratios":height_ratios, "width_ratios" : [3.,3.,3.,1.,1]}

def process_outlier(exclude,subset):
    for outlier in outliers:
        if exclude and outlier in subset.columns:
            subset = subset.loc[:,~subset.columns.str.contains('Tinaja Refed Liver 6')]
    return subset

fig,ax = plt.subplots(nrows=len(categories),ncols=5,figsize=(12, 8), gridspec_kw=gridspec_kw)
for i,category in zip(range(len(categories)),categories):
    for j,tissue in zip(range(3),['Brain', 'Muscle', 'Liver']):
        # get subset
        subset = process_outlier(args.exclude_outlier,DataFrame(astyanax_data.loc[
            ali.lmdata[cattype] == category,
            astyanax_data.columns.str.contains(tissue)]))

        pca_data = StandardScaler().fit_transform(subset.transpose())
        pca = PCA(n_components=2)
        pca_analysis = pca.fit(pca_data)
        tf_data = DataFrame(pca_analysis.transform(pca_data), columns=['C1','C2'], index=subset.columns)

        for color_theme,pop in zip(['Blues', 'Greens', "Reds"], ['Surface', 'Tinaja', 'Pachon']):
            cmap = matplotlib.cm.get_cmap(color_theme)
            for color_pos,feeding_state in zip([0.25, 0.5, 0.9], ['4d Starved', '30d Starved', 'Refed']):
                p = tf_data.iloc[subset.columns.str.contains(feeding_state) & subset.columns.str.contains(pop)]
                label = ', '.join((pop, feeding_state))
                ax[i,j].scatter(p['C1'],p['C2'], label=label, color=cmap(color_pos))
        handles, labels = ax[i,j].get_legend_handles_labels()

        if i == 0:
            ax[i,j].set_title(tissue)

for i,category in zip(range(len(categories)),categories):
    xpos = 0.1
    if category == 'Glycerophospholipids':
        category = 'Glycerophospho-\nlipids'
    ax[i,3].text(xpos,0.5,category,size=11,rotation=90.,ha='center',va='center',transform=ax[i,3].transAxes)
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

fig.legend(handles, labels, loc='center right')

plt.savefig(args.output)

