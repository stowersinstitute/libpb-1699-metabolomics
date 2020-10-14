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
parser.add_argument('--output', type=str, help="Output file.")
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


def tidy(data):
    output = []
    for pop in pops:
        for tissue in tissues:
            for condition in conditions:
                for i,row in data.loc[:,
                        data.columns.str.contains(pop) & data.columns.str.contains(tissue) & data.columns.str.contains(condition)].iterrows():
                    for k,u in enumerate(row,start=1):
                        wt = float(weight_matrix.loc[(weight_matrix['Population'] == pop) & (weight_matrix['Tissue'] == tissue) & (weight_matrix['Condition'] == condition),'Mass (mg)'].iloc[k-1])
                        output.append({'Name': i, 'Population': pop, 'Tissue': tissue, 'Condition': condition, 'Replicate': k, 'Weight': wt, 'Value': u})
    return DataFrame(output)

tidy_mtic = tidy(astyanax_data)
#print(tidy_mtic)
print(tidy_mtic.loc[tidy_mtic['Weight'] < 2.].groupby(['Population','Tissue']).median())
print(tidy_mtic.groupby(['Population','Tissue']).median())

# plot low weight distributions
low_wt_data = {}
for pop in ['Surface','Pachon']:
    for tissue in ['Liver']:
        for condition in conditions:
            low_wt_data[pop,tissue,condition] = astyanax_data.loc[:,astyanax_data.columns.str.contains(pop) & astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition) & list(weights_series < 2.)]

for i,pop in zip(range(len(pops)),pops):
    for tissue in tissues:
        if pop in ['Surface','Pachon'] and tissue in ['Liver']:
            d = concat((low_wt_data[pop,tissue,condition] for condition in conditions), axis=1)

unnormalized_data = concat((ame.row_table['KEGG'], read_csv(os.path.join(os.environ['JENNA_METABOLOMICS_PREFIX'],'unnormalized.csv'), skiprows=8).iloc[:, 8:-3]), axis=1).dropna()
unnormalized_data = unnormalized_data.set_index('KEGG')
unnormalized_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
unnormalized_data = unnormalized_data.loc[:,['pools' not in c for c in unnormalized_data.columns]]
unnormalized_data = unnormalized_data.apply(log10)
unnormalized_data = unnormalized_data.rename(ame.get_kegg_to_name_map(), axis=0)

tidy_unnorm = tidy(unnormalized_data)


w = weights_series.mean()
for k in range(len(unnormalized_data.columns)):
    unnormalized_data.iloc[:,k] += log10(w)-log10(weights_series.iloc[k])

tidy_wt_norm = tidy(unnormalized_data)

table = r'''
\begin{tabular}{ r | r | r}
\multicolumn{1}{c}{\textbf{Scheme}} & \multicolumn{1}{c}{\textbf{Tissue}} & \multicolumn{1}{c}{\textbf{Population}} & \multicolumn{1}{c}{\thead{Median \\ Low Wt. \\ Samples}} & \multicolumn{1}{c}{\textbf{Median \\ Rest}} \\ \hline
'''
def make_lines(input_data, scheme, tissue):
    lines = ''
    for pop in pops:
        subset = input_data.loc[(input_data['Population'] == pop) & (input_data['Tissue'] == tissue)].reset_index()
        lowwt = subset.loc[input_data['Weight'] < 2.].reset_index()
        lines += r'{scheme} & {tissue} & {pop} & ' + ' & '.join([d.groupby(['Population']).median()['Value'] for d in [subset,lowwt]]) + r'\\'
    return lines
table += make_lines(tidy_mtic, 'mTIC', 'Liver')
#print(tidy_mtic)
#print(tidy_unnorm)
#print(tidy_wt_norm)

table += '\n'+r'\end{tabular}'

print(table)
with open(args.output, 'w') as f:
    f.write(table)



