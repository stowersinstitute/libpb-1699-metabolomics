from argparse import ArgumentParser
import os

from numpy import log10, array, corrcoef, diag, cov, dot, log10, exp, power
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

from cavefinomics import AstyanaxMe, AstyanaxLi, MaMammalianMe, run_opls, plot_opls

parser = ArgumentParser(description="PCA vs mammals.")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument("--lipidmaps-json", type=str, help="Lipidmaps JSON.")
parser.add_argument("--output-dir", type=str, help="The output directory.")
args = parser.parse_args()

#ame = AstyanaxMe(
    #data_csv=args.astyanax,
    #kegg_compounds_file=args.compounds,
    #hmdb_file=args.hmdb,
    #)
ali = AstyanaxLi(
    lipids_normalized=args.lipids_normalized,
    lipidmaps_js=args.lipidmaps_json,
    )


#mammals = MaMammalianMe(annotation_csv=args.mammals_annotation, normalized_value_csv=args.mammals_normalized)

# get data and compute common metabolites
astyanax_data = ali.lmdata.iloc[:,ali.non_numeric_cols-1:]
astyanax_data.columns = (c[0] + ' ' + c[-2] + ' ' + ' '.join(c[1:-2]) + ' ' + c[-1] for c in (cc.split() for cc in astyanax_data.columns))
astyanax_data = astyanax_data.reindex(sorted(astyanax_data.columns), axis=1)
astyanax_data = astyanax_data.apply(log10)
#print(astyanax_data.columns)

#weights_series = ame.column_table['Mass (mg)']
#weights_series.index = [' '.join(u) for u in ame.treatment_descriptors]
#weights_series = weights_series.loc[['pools' not in c for c in weights_series.index]]
#astyanax_data = astyanax_data.apply(log10)

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = {'30d':'30d Starved', '4d':'4d Starved', 'Ref':'Refed'}
comparisons = {'PvS':('Pachon','Surface'),'TvS':('Tinaja','Surface'),'PvT':('Pachon','Tinaja')}
cattypes = {'Class':'class', 'Category':'category'}
outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def process_outlier(subset,comp):
    #if exclude_outlier and tissue=='Liver':
        #print(subset.columns)
        #if comp == 'TvS':
            #cols = list(subset.columns)
            #cols[11] = 'dropme'
            #subset.columns = cols
            #return subset.drop('dropme',axis=1)
        #elif comp == 'PvT':
            #cols = list(subset.columns)
            #cols[11] = 'dropme'
            #subset.columns = cols
            #return subset.drop('dropme',axis=1)
        #else:
            #return subset
    #else:
        #return subset
    for o in outliers:
        subset = subset.loc[:,~subset.columns.str.contains(o)]
    return subset

def process(subset, groups):
    normalized_data = DataFrame(StandardScaler().fit_transform(subset.transpose()), index=subset.columns, columns=subset.index)
    normalized_data.columns.name = 'InChIKey'
    # get pop name
    normalized_data = normalized_data.rename(lambda u: u.split()[0])
    #print(cattype,category,tissue,cond,comp)
    #if cattype == 'Class' and category == 'Docosanoids' and tissue == 'Brain' and cond == '30d' and comp == 'PvS':
        #print(normalized_data)
    if len(normalized_data.columns) > 1:
        pls, opls, Z, q2, dq2, p, acc, y_pred_pre, y_pred = run_opls(
          normalized_data,
          array([1 if groups[0] in u else -1 for u in subset.columns]),
          category_names=groups,
          population=None,
          tissue=None,
          n_components=1,
          n_p_val_iter=0,
          )
        nonortho_data = DataFrame(data=Z, index=normalized_data.index, columns=normalized_data.columns)
    else:
        nonortho_data = normalized_data
    nonortho_data.columns.name = 'InChIKey'
    return normalized_data, nonortho_data

for tissue in tissues:
    for cond,condition in conditions.items():
        for comp,groups in comparisons.items():
            for exclude_outlier,outlier_text in zip([False,True],['mit-Ausreißern','kein-Ausreißern']):
                for cattype in cattypes:
                    not_group = [c for c in pops if c not in groups][0]
                    #print(tissue,condition,not_group)
                    subset = DataFrame(astyanax_data.loc[
                        :,
                        astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition) &
                        ~astyanax_data.columns.str.contains(not_group)])
                    subset = process_outlier(subset,comp)
                    normalized_data, nonortho_data = process(subset, groups)

                    outdir = f"{args.output_dir}/opls/{outlier_text}/{tissue}/{cond}"
                    #https://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
                    # write opls data
                    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
                    fp = os.path.join(outdir,f'{comp}.csv')
                    nonortho_data.to_csv(fp)
                    #print(f'wrote {fp}')
                    # write z score-normalized data
                    outdir = f"{args.output_dir}/zscore/{outlier_text}/{tissue}/{cond}/"
                    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
                    fp = os.path.join(outdir,f'{comp}.csv')
                    normalized_data.to_csv(fp)
                    #print(f'wrote {fp}')

                    # categories
                    categories = list(sorted(set(ali.lmdata[cattype])))
                    for category in categories:
                        subset = DataFrame(astyanax_data.loc[
                            ali.lmdata[cattype] == category,
                            astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition) &
                            ~astyanax_data.columns.str.contains(not_group)])
                        subset = process_outlier(subset,comp)

                        normalized_data, nonortho_data = process(subset, groups)

                        c = category.replace(' ', '_').replace('\n','').replace('/','-')
                        outdir = f"{args.output_dir}/opls/{outlier_text}/{cattypes[cattype]}/{c}/{tissue}/{cond}"
                        #https://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
                        # write opls data
                        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
                        fp = os.path.join(outdir,f'{comp}.csv')
                        nonortho_data.to_csv(fp)
                        #print(f'wrote {fp}')
                        # write z score-normalized data
                        outdir = f"{args.output_dir}/zscore/{outlier_text}/{cattypes[cattype]}/{c}/{tissue}/{cond}/"
                        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
                        fp = os.path.join(outdir,f'{comp}.csv')
                        normalized_data.to_csv(fp)
                        #print(f'wrote {fp}')

