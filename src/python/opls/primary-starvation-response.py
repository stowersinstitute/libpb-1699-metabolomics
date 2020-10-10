from argparse import ArgumentParser
import os

from numpy import log, array, corrcoef, diag, cov, dot, log10, exp, power
from math import sqrt
from pandas import DataFrame, concat, Series
from sklearn.preprocessing import StandardScaler, Normalizer
import matplotlib.pyplot as plt
import matplotlib
from statsmodels.api import OLS
from statsmodels.formula.api import glm
import statsmodels.formula.api as smf
import statsmodels.api as sm
from seaborn import violinplot, boxplot, lmplot, regplot, scatterplot, lineplot, heatmap
from scipy.signal import correlate2d
from scipy.stats import ttest_ind
import scipy.stats
from pprint import pprint
import pathlib
from itertools import chain, repeat

from cavefinomics import AstyanaxMe, run_opls, plot_opls

parser = ArgumentParser(description="PCA vs mammals.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
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

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = ['30d Starved', '4d Starved', 'Refed']
comparisons = {'30vR':('30d Starved','Refed'),'4vR':('4d Starved','Refed'),'30v4':('30d Starved','4d Starved')}
categories = list(sorted(ame.compounds_by_category_from_dataset.keys()))
outliers = ['Pachon Muscle Refed 5','Tinaja Liver Refed 6','Pachon Liver 30d Starved 3']

for j,category in zip(range(len(categories)),categories):
    for i,tissue in zip(range(3),['Brain', 'Muscle', 'Liver']):
        for pop in ['Pachon', 'Tinaja', 'Surface']:
            for comp,groups in comparisons.items():
                for exclude_outlier,outlier_text in zip([False,True],['with-outliers','without-outliers']):
                    not_group = [c for c in conditions if c not in groups][0]
                    subset = DataFrame(astyanax_data.iloc[
                        astyanax_data.index.isin(ame.compounds_by_category_from_dataset[category]),
                        astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(pop) &
                        ~astyanax_data.columns.str.contains(not_group)])
                    subset.columns = [' '.join((c,str(n))) for c,n in zip(subset.columns,chain.from_iterable(repeat(range(1,7),2)))]
                    if exclude_outlier:
                        for outlier in outliers:
                            if outlier in subset.columns:
                                subset = subset.drop(outlier,axis=1)
                    subset.columns = [' '.join(c.split()[2:-1]) for c in subset.columns]
                    kegg = list(subset.index)
                    subset = subset.rename(ame.get_kegg_to_name_map(), axis=0)

                    normalized_data = DataFrame(StandardScaler().fit_transform(subset.transpose()), index=subset.columns, columns=subset.index)
                    pls, opls, Z, q2, dq2, p, acc, y_pred_pre, y_pred = run_opls(
                      normalized_data,
                      array([1 if groups[0] in u else -1 for u in subset.columns]),
                      category_names=groups,
                      population=pop,
                      tissue=tissue,
                      n_components=1,
                      n_p_val_iter=0,
                      )
                    nonortho_data = DataFrame(data=Z, index=normalized_data.index, columns=normalized_data.columns)
                    nonortho_data.columns.name = 'Name'
                    s = DataFrame(kegg, index=nonortho_data.columns, columns=['KEGG']).transpose()
                    nonortho_data = concat((s,nonortho_data),axis=0)

                    c = category.replace(' ', '_').replace('\n','').replace('/','-')
                    outdir = f"{args.output_dir}/opls/{outlier_text}/{c}/{tissue}/{pop}"
                    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
                    nonortho_data.to_csv(os.path.join(outdir,f'{comp}.csv'))

