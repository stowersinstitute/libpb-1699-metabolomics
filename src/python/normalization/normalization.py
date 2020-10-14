from argparse import ArgumentParser
import os

from numpy import array, log10, c_
from pandas import DataFrame, Series, concat, read_csv
from sklearn.preprocessing import StandardScaler, Normalizer
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import matplotlib
from seaborn import heatmap, distplot
from sklearn.decomposition import PCA as PCA
from functools import reduce
from scipy.stats import chisquare, fisher_exact
import math
from pprint import pprint
from collections import defaultdict
import json
import seaborn as sns
matplotlib.style.use('seaborn')

from cavefinomics import AstyanaxMe

#https://stackoverflow.com/questions/34706845/change-xticklabels-fontsize-of-seaborn-heatmap
sns.set(font_scale=1.4)

parser = ArgumentParser(description="Plot distribution of peak intensities for various samples.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--unnormalized", type=str, help="Unnormalized metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument('--output-dir', type=str, help="Output directory.")
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

# calculate weight matrix

weights = {}
weight_matrix = []
for pop in pops:
    for tissue in tissues:
        for condition in conditions:
            weights[pop,tissue,condition] = weights_series.loc[
                    astyanax_data.columns.str.contains(pop) & astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition)]
            for i,v in weights[pop,tissue,condition].iteritems():
                weight_matrix.append({'Population':pop,'Tissue':tissue,'Condition':condition,'Mass (mg)':v})

weight_matrix = DataFrame(weight_matrix)


all_pachon_liver = weight_matrix.loc[weight_matrix['Population'].str.contains('Pachon') & weight_matrix['Tissue'].str.contains('Liver')]
pachon_liver_low_wt = all_pachon_liver.loc[weight_matrix['Mass (mg)'] < 2.]
#https://stackoverflow.com/questions/24644656/how-to-print-pandas-dataframe-without-index
all_surface_liver = weight_matrix.loc[weight_matrix['Population'].str.contains('Surface') & weight_matrix['Tissue'].str.contains('Liver')]
surface_liver_low_wt = all_surface_liver.loc[weight_matrix['Mass (mg)'] < 2.]

subsets = {}
for pop in pops:
    for tissue in tissues:
        for condition in conditions:
            subsets[pop,tissue,condition] = astyanax_data.loc[:,
                    astyanax_data.columns.str.contains(pop) & astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition)]

gridspec_kw = {"height_ratios":[1.]*len(pops), "width_ratios" : [3.,3.,3.,3.,1]}
fig,ax = plt.subplots(nrows=len(pops),ncols=len(tissues)+2,sharex=False,sharey=False,gridspec_kw=gridspec_kw,figsize=(12.,8.))
fig.suptitle('$\log_{10}$ mTIC data')

for i,pop,color in zip(range(len(pops)),pops,['firebrick','goldenrod','dodgerblue']):
    for j,tissue in enumerate(tissues):
        d = concat((subsets[pop,tissue,condition] for condition in conditions), axis=1)
        distplot(d.values.flatten(),color=color,ax=ax[i,j])
        ax[i,j].set_xlim((0.,7.5))
        ax[i,j].set_ylim((0.,.5))
        if i == 0:
            ax[i,j].set_title(tissue)
        if j>0:
            ax[i,j].yaxis.set_ticklabels([])
        if i+1<3:
            #https://stackoverflow.com/questions/2176424/hiding-axis-text-in-matplotlib-plots
            ax[i,j].xaxis.set_ticklabels([])
    ax[i,4].text(0.5,0.5,pop,size=16,rotation=90.,ha='center',va='center',transform=ax[i,4].transAxes)
    # hide graphics
    ax[i,4].set_yticks([], minor=[])
    ax[i,4].set_xticks([], minor=[])
    ax[i,4].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        ax[i,4].spines[s].set_visible(False)

# plot low weight distributions
low_wt_data = {}
for pop in ['Surface','Pachon']:
    for tissue in ['Liver']:
        for condition in conditions:
            low_wt_data[pop,tissue,condition] = astyanax_data.loc[:,astyanax_data.columns.str.contains(pop) & astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition) & list(weights_series < 2.)]

for i,pop,color in zip(range(len(pops)),pops,['firebrick','goldenrod','dodgerblue']):
    for tissue in tissues:
        if pop in ['Surface','Pachon'] and tissue in ['Liver']:
            d = concat((low_wt_data[pop,tissue,condition] for condition in conditions), axis=1)
            distplot(d.values.flatten(),color=color,ax=ax[i,3])
            ax[i,3].set_xlim((0.,7.5))
            ax[i,3].set_ylim((0.,.5))
        ax[i,3].yaxis.set_ticklabels([])
        if i+1<3:
            ax[i,3].xaxis.set_ticklabels([])

# remove axis labels
for i in range(len(pops)):
    for j in [1,2,3,4]:
        ax[i,j].set_ylabel('')

ax[0,3].set_title('Low Wt Liver')
ax[1,3].text(0.5,0.5,'N/A',size=16,ha='center',va='center',transform=ax[1,3].transAxes)

plt.savefig(os.path.join(args.output_dir,'density-mtic-normalized.pdf'))

unnormalized_data = concat((ame.row_table['KEGG'], read_csv(args.unnormalized, skiprows=8).iloc[:, 8:-3]), axis=1).dropna()
unnormalized_data = unnormalized_data.set_index('KEGG')
unnormalized_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
unnormalized_data = unnormalized_data.loc[:,['pools' not in c for c in unnormalized_data.columns]]
unnormalized_data = unnormalized_data.apply(log10)


unnormalized_data = unnormalized_data.rename(ame.get_kegg_to_name_map(), axis=0)


# **
# redo plot with tissue unnormalized data
# **

subsets = {}
for pop in pops:
    for tissue in tissues:
        for condition in conditions:
            subsets[pop,tissue,condition] = unnormalized_data.loc[:,
                    unnormalized_data.columns.str.contains(pop) & unnormalized_data.columns.str.contains(tissue) & unnormalized_data.columns.str.contains(condition)]

gridspec_kw = {"height_ratios":[1.]*len(pops), "width_ratios" : [3.,3.,3.,3.,1]}
fig,ax = plt.subplots(nrows=len(pops),ncols=len(tissues)+2,sharex=False,sharey=False,gridspec_kw=gridspec_kw,figsize=(12.,8.))
fig.suptitle('$\log_{10}$ Unnormalized Data')

for i,pop,color in zip(range(len(pops)),pops,['firebrick','goldenrod','dodgerblue']):
    for j,tissue in enumerate(tissues):
        d = concat((subsets[pop,tissue,condition] for condition in conditions), axis=1)
        distplot(d.values.flatten(),color=color,ax=ax[i,j])
        ax[i,j].set_xlim((0.,7.5))
        ax[i,j].set_ylim((0.,.5))
        if i == 0:
            ax[i,j].set_title(tissue)
        if j>0:
            ax[i,j].yaxis.set_ticklabels([])
        if i+1<3:
            #https://stackoverflow.com/questions/2176424/hiding-axis-text-in-matplotlib-plots
            ax[i,j].xaxis.set_ticklabels([])
    ax[i,4].text(0.5,0.5,pop,size=16,rotation=90.,ha='center',va='center',transform=ax[i,4].transAxes)
    # hide graphics
    ax[i,4].set_yticks([], minor=[])
    ax[i,4].set_xticks([], minor=[])
    ax[i,4].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        ax[i,4].spines[s].set_visible(False)

# plot low weight distributions
low_wt_data = {}
for pop in ['Surface','Pachon']:
    for tissue in ['Liver']:
        for condition in conditions:
            low_wt_data[pop,tissue,condition] = unnormalized_data.loc[:,unnormalized_data.columns.str.contains(pop) & unnormalized_data.columns.str.contains(tissue) & unnormalized_data.columns.str.contains(condition) & list(weights_series < 2.)]

for i,pop,color in zip(range(len(pops)),pops,['firebrick','goldenrod','dodgerblue']):
    for tissue in tissues:
        if pop in ['Surface','Pachon'] and tissue in ['Liver']:
            d = concat((low_wt_data[pop,tissue,condition] for condition in conditions), axis=1)
            distplot(d.values.flatten(),color=color,ax=ax[i,3])
            ax[i,3].set_xlim((0.,7.5))
            ax[i,3].set_ylim((0.,.5))
        ax[i,3].yaxis.set_ticklabels([])
        if i+1<3:
            ax[i,3].xaxis.set_ticklabels([])

# remove axis labels
for i in range(len(pops)):
    for j in [1,2,3,4]:
        ax[i,j].set_ylabel('')

ax[0,3].set_title('Low Wt Liver')
ax[1,3].text(0.5,0.5,'N/A',size=16,ha='center',va='center',transform=ax[1,3].transAxes)

plt.savefig(os.path.join(args.output_dir,'density-unnormalized.pdf'))

# **
# redo plot with tissue weight-normalized data
# **


w = weights_series.mean()
for k in range(len(unnormalized_data.columns)):
    unnormalized_data.iloc[:,k] += log10(w)-log10(weights_series.iloc[k])

subsets = {}
for pop in pops:
    for tissue in tissues:
        for condition in conditions:
            subsets[pop,tissue,condition] = unnormalized_data.loc[:,
                    unnormalized_data.columns.str.contains(pop) & unnormalized_data.columns.str.contains(tissue) & unnormalized_data.columns.str.contains(condition)]

gridspec_kw = {"height_ratios":[1.]*len(pops), "width_ratios" : [3.,3.,3.,3.,1]}
fig,ax = plt.subplots(nrows=len(pops),ncols=len(tissues)+2,sharex=False,sharey=False,gridspec_kw=gridspec_kw,figsize=(12.,8.))
fig.suptitle('$\log_{10}$ Sample Weight-normalized Data')

for i,pop,color in zip(range(len(pops)),pops,['firebrick','goldenrod','dodgerblue']):
    for j,tissue in enumerate(tissues):
        d = concat((subsets[pop,tissue,condition] for condition in conditions), axis=1)
        distplot(d.values.flatten(),color=color,ax=ax[i,j])
        ax[i,j].set_xlim((0.,7.5))
        ax[i,j].set_ylim((0.,.5))
        if i == 0:
            ax[i,j].set_title(tissue)
        if j>0:
            ax[i,j].yaxis.set_ticklabels([])
        if i+1<3:
            #https://stackoverflow.com/questions/2176424/hiding-axis-text-in-matplotlib-plots
            ax[i,j].xaxis.set_ticklabels([])
    ax[i,4].text(0.5,0.5,pop,size=16,rotation=90.,ha='center',va='center',transform=ax[i,4].transAxes)
    # hide graphics
    ax[i,4].set_yticks([], minor=[])
    ax[i,4].set_xticks([], minor=[])
    ax[i,4].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        ax[i,4].spines[s].set_visible(False)

# plot low weight distributions
low_wt_data = {}
for pop in ['Surface','Pachon']:
    for tissue in ['Liver']:
        for condition in conditions:
            low_wt_data[pop,tissue,condition] = unnormalized_data.loc[:,unnormalized_data.columns.str.contains(pop) & unnormalized_data.columns.str.contains(tissue) & unnormalized_data.columns.str.contains(condition) & list(weights_series < 2.)]

for i,pop,color in zip(range(len(pops)),pops,['firebrick','goldenrod','dodgerblue']):
    for tissue in tissues:
        if pop in ['Surface','Pachon'] and tissue in ['Liver']:
            d = concat((low_wt_data[pop,tissue,condition] for condition in conditions), axis=1)
            distplot(d.values.flatten(),color=color,ax=ax[i,3])
            ax[i,3].set_xlim((0.,7.5))
            ax[i,3].set_ylim((0.,.5))
        ax[i,3].yaxis.set_ticklabels([])
        if i+1<3:
            ax[i,3].xaxis.set_ticklabels([])

# remove axis labels
for i in range(len(pops)):
    for j in [1,2,3,4]:
        ax[i,j].set_ylabel('')

ax[0,3].set_title('Low Wt Liver')
ax[1,3].text(0.5,0.5,'N/A',size=16,ha='center',va='center',transform=ax[1,3].transAxes)

plt.savefig(os.path.join(args.output_dir,'density-weight-normalized.pdf'))
