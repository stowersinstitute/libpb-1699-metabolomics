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

from cavefinomics import AstyanaxMe, MaMammalianMe

parser = ArgumentParser(description="PCA for primary.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--mammals-annotation", type=str, help="Annotation for mammals.")
parser.add_argument("--mammals-normalized", type=str, help="Normalized mammalian metabolite file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--groups", type=str, help="Groups to compare (30d Starved, 4d Starved, Refed).")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the single Tinaja outlier?")
parser.add_argument("--output", type=str, help="Output primary PCA.")
args = parser.parse_args()

short_groups = ['4d Starved', 'Refed']
long_groups = ['30d Starved', 'Refed']

ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    hmdb_file=args.hmdb,
    )

mammals = MaMammalianMe(annotation_csv=args.mammals_annotation, normalized_value_csv=args.mammals_normalized)

# get data and compute common metabolites
astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
astyanax_data = astyanax_data.apply(log10)


categories = list(sorted(ame.compounds_by_category_from_dataset.keys()))
outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def process_outlier(exclude,subset):
    for outlier in outliers:
        if exclude and outlier in subset.columns:
            subset = subset.loc[:,~subset.columns.str.contains('Tinaja Refed Liver 6')]
    return subset

if True:
    height_ratios = [1.]*len(categories)
else:
    height_ratios = [float(len(ame.compounds_by_category_from_dataset[c])) for c in categories]
gridspec_kw = {"height_ratios":height_ratios, "width_ratios" : [3.,3.,3.,1.,1]}

f,axs = plt.subplots(nrows=len(categories),ncols=5,figsize=(12, 8), gridspec_kw=gridspec_kw)
for j,category in zip(range(len(categories)),categories):
    for i,tissue in zip(range(3),['Brain', 'Muscle', 'Liver']):
        # get 4d starved, 30d starved, and refed groups
        subset = process_outlier(args.exclude_outlier,DataFrame(astyanax_data.iloc[
            astyanax_data.index.isin(ame.compounds_by_category_from_dataset[category]),
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
                axs[j,i].scatter(p['C1'],p['C2'], label=label, color=cmap(color_pos))
        handles, labels = axs[j,i].get_legend_handles_labels()

        if j == 0:
            axs[j,i].set_title(tissue)

for j,category in zip(range(len(categories)),categories):
    xpos = 0.1
    axs[j,3].text(xpos,0.5,category,size=11,rotation=90.,ha='center',va='center',transform=axs[j,3].transAxes)
    # hide graphics
    axs[j,3].set_yticks([], minor=[])
    axs[j,3].set_xticks([], minor=[])
    axs[j,3].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        axs[j,3].spines[s].set_visible(False)
    # hide graphics
    axs[j,4].set_yticks([], minor=[])
    axs[j,4].set_xticks([], minor=[])
    axs[j,4].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        axs[j,4].spines[s].set_visible(False)


f.legend(handles, labels, loc='center right')

plt.savefig(args.output)

