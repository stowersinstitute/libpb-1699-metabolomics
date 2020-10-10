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

from cavefinomics import AstyanaxMe

parser = ArgumentParser(description="PCA for primary.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--groups", type=str, help="Groups to compare (30d Starved, 4d Starved, Refed).")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the single Tinaja outlier?")
parser.add_argument("--output", type=str, help="Output primary PCA.")
args = parser.parse_args()

if not args.exclude_outlier:
    raise RuntimeError('should exclude outlier')

short_groups = ['4d Starved', 'Refed']
long_groups = ['30d Starved', 'Refed']

ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    sample_sheet_path=args.sample_sheet,
    hmdb_file=args.hmdb,
    )


# get data and compute common metabolites
astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
astyanax_data = astyanax_data.apply(log10)


categories = list(sorted(ame.compounds_by_category_from_dataset.keys()))
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

def add_replicate_nums(subset):
    subset.columns = (f'{c} {n}' for n,c in zip(list(range(1,7))*9,subset.columns))
    return subset

def process_outlier(exclude,subset):
    for outlier in outliers:
        if exclude and outlier in subset.columns:
            subset = subset.loc[:,~subset.columns.str.contains(outlier)]
    return subset


remap_cond = {
  '4d Starved': '4d Fasted',
  '30d Starved': '30d Fasted',
  }

fig,ax = plt.subplots(figsize=(8, 8))
for i,category in zip(range(len(categories)),categories):
    for j,tissue in zip(range(3),['Brain', 'Muscle', 'Liver']):
        # get 4d starved, 30d starved, and refed groups
        subset = process_outlier(args.exclude_outlier,add_replicate_nums(DataFrame(astyanax_data.iloc[
            astyanax_data.index.isin(ame.compounds_by_category_from_dataset[category]),
            astyanax_data.columns.str.contains(tissue)])))

        pca_data = StandardScaler().fit_transform(subset.transpose())
        pca = PCA(n_components=2)
        pca_analysis = pca.fit(pca_data)
        tf_data = DataFrame(pca_analysis.transform(pca_data), columns=['C1','C2'], index=subset.columns)

        ax.cla()
        for pop in ['Surface', 'Tinaja', 'Pachon']:
            for color,feeding_state in zip([adjust_lightness(colors[pop],c) for c in [0.5, 1.0, 1.5]], ['4d Starved', '30d Starved', 'Refed']):
                p = tf_data.iloc[subset.columns.str.contains(feeding_state) & subset.columns.str.contains(pop)]
                label = ', '.join((pop, remap_cond[feeding_state] if feeding_state in remap_cond else feeding_state))
                ax.scatter(p['C1'],p['C2'], label=label, color=color)
        handles, labels = ax.get_legend_handles_labels()

        if i == 0:
            ax.set_title(tissue)

fig,ax = plt.subplots(figsize=(2.5, 2.5))

# hide graphics
ax.set_yticks([], minor=[])
ax.set_xticks([], minor=[])
ax.patch.set_visible(False)
for s in ["top", "bottom", "left", "right"]:
    ax.spines[s].set_visible(False)

fig.legend(handles, labels, loc='center')

plt.savefig(args.output)

