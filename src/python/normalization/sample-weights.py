from argparse import ArgumentParser
import os

from numpy import array, log10, c_
import numpy as np
from pandas import DataFrame, Series, concat
from sklearn.preprocessing import StandardScaler, Normalizer
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import matplotlib
from seaborn import heatmap, distplot, stripplot, pointplot
import seaborn as sns
from sklearn.decomposition import PCA as PCA
from functools import reduce
from scipy.stats import chisquare, fisher_exact
import math
from pprint import pprint
from collections import defaultdict
import json

from cavefinomics import AstyanaxMe

#https://stackoverflow.com/questions/34706845/change-xticklabels-fontsize-of-seaborn-heatmap
sns.set(font_scale=1.2)
matplotlib.style.use('default')

parser = ArgumentParser(description="Plot tissue sample weights.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--output", type=str, help="Output.")
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
astyanax_data = astyanax_data.apply(log10)


categories = list(sorted(ame.compounds_by_category_from_dataset.keys()))

astyanax_data = astyanax_data.rename(ame.get_kegg_to_name_map(), axis=0)

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = ['30d Starved', '4d Starved', 'Refed']

# stabilize random numbers
#https://stackoverflow.com/questions/61944815/how-to-set-seed-for-jitter-in-seaborn-stripplot
np.random.seed(101)

weights = {}
weight_matrix = []
for pop in pops:
    for tissue in tissues:
        for condition in conditions:
            weights[pop,tissue,condition] = weights_series.loc[
                    astyanax_data.columns.str.contains(pop) & astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition)]
            #print(weights[pop,tissue,condition])
            for i,v in weights[pop,tissue,condition].iteritems():
                weight_matrix.append({'Population':pop,'Tissue':tissue,'Condition':condition,'Mass (mg)':v})

weight_matrix = DataFrame(weight_matrix)

fig,ax = plt.subplots()
sns.despine(bottom=True, left=True)


all_pachon_liver = weight_matrix.loc[weight_matrix['Population'].str.contains('Pachon') & weight_matrix['Tissue'].str.contains('Liver')]
pachon_liver_low_wt = all_pachon_liver.loc[weight_matrix['Mass (mg)'] < 2.]
#https://stackoverflow.com/questions/24644656/how-to-print-pandas-dataframe-without-index
print(pachon_liver_low_wt.to_string(index=False))
all_surface_liver = weight_matrix.loc[weight_matrix['Population'].str.contains('Surface') & weight_matrix['Tissue'].str.contains('Liver')]
surface_liver_low_wt = all_surface_liver.loc[weight_matrix['Mass (mg)'] < 2.]
print(surface_liver_low_wt.to_string(index=False))

stripplot(x="Mass (mg)", y="Tissue", hue="Population", palette=['firebrick','goldenrod','dodgerblue'],
              data=weight_matrix, dodge=True, alpha=.5, zorder=1)

ax.annotate(
    f'{len(pachon_liver_low_wt.index)} / {len(all_pachon_liver.index)} samples < 2 mg\nin Pachon liver',
    xy=(0.165, 0.25),
    xycoords="axes fraction",
    xytext=(0.45, 0.45),
    textcoords="axes fraction",
    arrowprops=dict(facecolor="black", shrink=0.05, linewidth=0.0, alpha=0.75),
    horizontalalignment="center",
    verticalalignment="top",
    bbox=dict(fc="white", alpha=0.5, boxstyle="round"),
)

ax.annotate(
    f'{len(surface_liver_low_wt.index)} / {len(all_surface_liver.index)} samples < 2 mg\nin surface liver',
    xy=(0.185, 0.07),
    xycoords="axes fraction",
    xytext=(0.65, 0.25),
    textcoords="axes fraction",
    arrowprops=dict(facecolor="black", shrink=0.05, linewidth=0.0, alpha=0.75),
    horizontalalignment="center",
    verticalalignment="top",
    bbox=dict(fc="white", alpha=0.5, boxstyle="round"),
)


plt.savefig(args.output)
