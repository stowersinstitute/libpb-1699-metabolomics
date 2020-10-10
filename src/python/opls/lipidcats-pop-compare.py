from argparse import ArgumentParser
from pandas import read_csv
from cavefinomics import AstyanaxLi
from seaborn import heatmap, distplot
from seaborn import catplot
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import matplotlib
from numpy import log10, isfinite, any, all, quantile, nan, isnan
from pandas import DataFrame, concat, melt, option_context
import pathlib
import os
import json
from cavefinomics import run_opls
from sklearn.preprocessing import StandardScaler, Normalizer
from numpy import log10, mean, isnan, quantile, nan, array

parser = ArgumentParser(description="PCA vs mammals.")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument("--lipidmaps-json", type=str, help="Lipidmaps JSON.")
parser.add_argument("--lipidmaps-fa", type=str, help="Lipidmaps fa.")
parser.add_argument("--output-dir", type=str, help="The output directory.")
args = parser.parse_args()

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
polarities = ['positive','negative']
conditions = {'30d':'30d Starved', '4d':'4d Starved', 'Ref':'Refed'}

with open(args.lipidmaps_fa) as f:
  lmfa = json.load(f)

ali = AstyanaxLi(
  lipids_normalized=args.lipids_normalized,
  lipidmaps_js=args.lipidmaps_json,
  )

def get_lipidmap_id(inchikey):
    if inchikey in ali.lipidmaps_inchikey:
        return ali.lipidmaps_inchikey[inchikey]
    elif '-'.join(inchikey.split('-')[:2]) in ali.lipidmaps_inchikey2:
        return ali.lipidmaps_inchikey2['-'.join(inchikey.split('-')[:2])]
    elif inchikey.split('-')[0] in ali.lipidmaps_inchikey1:
        return ali.lipidmaps_inchikey1[inchikey.split('-')[0]]
    return nan


outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def process_outlier(subset):
    for o in outliers:
        subset.loc[:,subset.columns.str.contains(o)] = nan
    return subset

for include_outliers,outlier_label in zip([True,False],['with-outliers','without-outliers']):
    cats = []
    classes = []
    def fix_cols(df,tissue):
        d = df.copy()
        def j(c):
            return ' '.join(c)
        d.columns = (f'{c[0]} {tissue} {j(c[1:])} {n}' for c,n in zip((cc.split() for cc in d.columns),list(range(1,7))*27))
        if include_outliers:
            d = process_outlier(d)
        d = d.reindex(sorted(d.columns), axis=1)
        return d

    for tissue in tissues:
        for polarity in polarities:
            lmd = DataFrame(ali.normalized[tissue,polarity].apply(lambda u: get_lipidmap_id(u.name),axis=1), columns=['LMID'])
            lmd['CATEGORY'] = lmd.apply(lambda u: ali.lipidmaps[u[0]]['CATEGORY'] if isinstance(u[0],str) else nan,axis=1)
            def assignclass(u):
                if isinstance(u[0],str):
                    lmid = u[0]
                    mc = ali.lipidmaps[lmid]['MAIN_CLASS']
                    return mc
                else:
                    return nan
            lmd['MAIN_CLASS'] = lmd.apply(assignclass,axis=1)
            d = ali.normalized[tissue,polarity]
            d = fix_cols(d,tissue)
            d['CATEGORY'] = lmd['CATEGORY']
            #https://github.com/pandas-dev/pandas/issues/29481#issuecomment-652647458
            sum_by_category = d.groupby('CATEGORY').apply(lambda x: x.sum(skipna=False))

            d = ali.normalized[tissue,polarity]
            d = fix_cols(d,tissue)
            d['MAIN_CLASS'] = lmd['MAIN_CLASS']
            #https://github.com/pandas-dev/pandas/issues/29481#issuecomment-652647458
            sum_by_class = d.groupby('MAIN_CLASS').apply(lambda x: x.sum(skipna=False))

            for i,condition in zip(range(3),conditions.values()):
                for j,pop in zip(range(1,4),pops):
                    category_subset = sum_by_category.loc[:,sum_by_category.columns.str.contains(condition) & sum_by_category.columns.str.contains(pop)]
                    classes_subset = sum_by_class.loc[:,sum_by_class.columns.str.contains(condition) & sum_by_class.columns.str.contains(pop)]

                    def process_df(d):
                        d = d.reset_index()
                        d.columns = ['Category'] + ['Intensity']*(len(d.columns)-1)
                        d.insert(0,'Tissue',[tissue]*len(d.index))
                        d.insert(1,'Polarity',[polarity]*len(d.index))
                        d.insert(2,'Population',[pop]*len(d.index))
                        d.insert(3,'Condition',[condition]*len(d.index))
                        d['Category'] = d['Category'].apply(lambda u: u.split('[')[0].strip())
                        return d

                    category_subset,classes_subset = (process_df(d) for d in [category_subset,classes_subset])

                    #print(category_subset)
                    cats.append(category_subset)
                    classes.append(classes_subset)

    cats = concat(cats,axis=0)
    classes = concat(classes,axis=0)
    cats = cats.groupby(['Tissue','Category','Population','Condition']).sum().reset_index()
    classes = classes.groupby(['Tissue','Category','Population','Condition']).sum().reset_index()

    cattypes = ['Categories','Classes']
    comparisons = {'PvS':['Pachon','Surface'],'TvS':['Tinaja','Surface'],'PvT':['Pachon','Tinaja']}

    dq2s = []
    for tissue in tissues:
        for condlabel,condition in conditions.items():
            for d,c in zip([cats,classes],cattypes):
                subset = {}
                for pop in pops:
                    df = d.loc[(d['Tissue'] == tissue) & (d['Population'] == pop) & (d['Condition'] == condition)]
                    df = df.set_index('Category')
                    df = df.iloc[:,3:]
                    df.columns = [pop]*len(df.columns)
                    subset[pop] = df

                for comp,groups in comparisons.items():
                    df = concat((subset[groups[0]].transpose(),subset[groups[1]].transpose()),axis=0)
                    df = df.apply(log10)
                    with option_context('mode.use_inf_as_null', True):
                        df = df.dropna()
                    normalized_data = DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns)
                    pls, opls, Z, q2, dq2, p, acc, y_pred_pre, y_pred = run_opls(
                          normalized_data,
                          array([1 if groups[0] in u else -1 for u in df.index]),
                          category_names=[groups[0], groups[1]],
                          population=pop,
                          tissue='any',
                          n_components=1,
                          n_p_val_iter=0,
                          )
                    dq2s.append({'Tissue':tissue,'Condition':condition,'Comparison':comp,'CatType':c,'DQ2':dq2})
                    pathname = f'{args.output_dir}/opls/{outlier_label}/{tissue}/{condlabel}/{c}'
                    pathlib.Path(pathname).mkdir(parents=True, exist_ok=True)
                    normalized_data.to_csv(os.path.join(pathname,f'{comp}.csv'),header=True)
    dq2s = DataFrame(dq2s)

