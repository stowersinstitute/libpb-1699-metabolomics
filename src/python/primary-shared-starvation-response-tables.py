from pandas import read_csv
from itertools import chain, zip_longest
from numpy import log2
import os
from argparse import ArgumentParser

parser = ArgumentParser(description="Adaptive response table.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--out-dir", type=str, help="Out dir")
args = parser.parse_args()

if args.exclude_outlier:
    outlier = 'kein-Ausreißern'
else:
    outlier = 'mit-Ausreißern'

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = ['30d Starved', '4d Starved', 'Refed']
comparisons = {'30vR':('30d Starved','Refed'),'4vR':('4d Starved','Refed'),'30v4':('30d Starved','4d Starved')}
categories = {"Aminoacids":'Amino acids',"Carbohydrates_-CCM": 'Carbohydrates / CCM',"Fattyacids":'Fatty acids',"Misc._-_sec.metabolites":'Misc. / sec. metabolites',"Nucleotides":'Nucleotides'}

sigmetabs = {}
for catname,category in categories.items():
    for comp,groups in comparisons.items():
        def get_sig_metabs(tissue):
            d = read_csv(f'out/work/primary/glm/singlefactor/{outlier}/{catname}/{tissue}/CvS/{comp}.csv')
            d = d.loc[d['Pr(>|z|)']<0.05]
            d = d.rename({'Pr(>|z|)':'p'},axis=1)
            cols = list(d.columns)
            cols[0] = 'Compound'
            d.columns = cols
            def fix_compound(c):
                c = c.replace('X2.hydroxyglutaric.acid','2-hydroxyglutaric acid')
                c = c.replace('X4.hydroxybenzoate','4-hydroxybenzoate')
                c = c.replace('X2.ketoadipic.acid','2-ketoadipic acid')
                c = c.replace('UDP.N.acetylglucosamine','UDP-N-acetylglucosamine')
                c = c.replace('.',' ')
                return c
            d['Compound'] = d['Compound'].apply(lambda c: fix_compound(c))
            d['Down'] = d['Estimate'] < 0.
            d = d.sort_values(['Down','p'])
            return d
        for tissue in tissues:
            sigmetabs[tissue,comp,category] = get_sig_metabs(tissue)

merged_data = read_csv('out/work/primary/merged-mtic.csv',index_col=0)
if args.exclude_outlier:
    merged_data = merged_data.loc[merged_data['Outlier'] == False]
merged_data = merged_data.drop(['KEGG','HMDB','ChEBI','Outlier'],axis=1)

catsub = {category:category for category in categories.values()}
catsub['Amino acids'] = r'\shortstack[l]{Amino\\ acids}'
catsub['Carbohydrates / CCM'] = r'\shortstack[l]{Carbo-\\hydrates /\\ CCM}'
catsub['Misc. / sec. metabolites'] = 'Misc.'
catlines = {}
for comp,groups in comparisons.items():
    for catname,category in categories.items():
        m = {}
        for tissue in tissues:
            group_data = [
              merged_data.loc[
                (merged_data['Tissue'] == tissue) &
                (merged_data['Category'] == category) &
                (merged_data['Condition'] == g)].drop(
                  'Category',axis=1) for g in groups]
            group_data = [g.groupby(['Name','Population','Tissue','Condition']).mean().drop(['Replicate'],axis=1).reset_index() for g in group_data]
            group_merged = group_data[0].merge(group_data[1],on=['Name','Population','Tissue'])
            group_merged['log2fc'] = log2(group_merged['Raw_mTIC_x']/group_merged['Raw_mTIC_y'])
            avg_log2fc = {cpd.strip():
                            0.5*float(group_merged.loc[(group_merged['Name'] == cpd) & (group_merged['Population'] == 'Pachon'),'log2fc']) +
                            0.5*float(group_merged.loc[(group_merged['Name'] == cpd) & (group_merged['Population'] == 'Tinaja'),'log2fc']) -float(group_merged.loc[(group_merged['Name'] == cpd) & (group_merged['Population'] == 'Surface'),'log2fc'])
                            for cpd in sorted(list(group_merged['Name']))}
            m[tissue] = list(
              zip(
                sigmetabs[tissue,comp,category]['Compound'],
                (f'{v:.2f}' for v in (avg_log2fc[cpd.strip()] for cpd in sigmetabs[tissue,comp,category]['Compound'])),
                )
              )
        lines = []
        k = 1
        l = max([len(m[t]) for t in tissues])
        for brain_vals,muscle_vals,liver_vals in zip_longest(m['Brain'],m['Muscle'],m['Liver'],fillvalue=('','')):
            if k == 1:
                firstrow = f'\\multirow{{{l}}}{{*}}{{{catsub[category]}}}'
            else:
                firstrow = ''
            k += 1
            lines.append(' & '.join([firstrow]+list(chain(brain_vals,muscle_vals,liver_vals))))
        catlines[comp,category] = '\n'.join(l+r' \\' for l in lines)

captions = {
  '30vR': 'Adaptive Response in 30-day Fasting',
  '4vR': 'Adaptive Response in 4-day Fasting',
  '30v4': 'Adaptive Response in 30- vs 4-day Fasting',
}

descs = {
  '30vR': r"Differentially significant metabolites in 30-day fasted states which are similar in both cave populations. This table shows p--values associated with a logistic regression model using OPLS--filtered z--scores for 30-day tasted and refed fish as input. To identify metabolites conserved between both cave populations, we implemented the test $0.5 \cdot (\mathrm{PS}-\mathrm{PR}+\mathrm{TS}-\mathrm{TR}) - (\mathrm{SS} - \mathrm{SR})$ where $\mathrm{PS}$ refers to Pach\'{o}n 30-day fasted, $\mathrm{PR}$ refers to Pach\'{o}n refed, etc. Tables \ref{table:significant-4vR} and \ref{table:significant-30v4} show conserved metabolites in cave populations versus surface for 4-day fasted versus refed and 30-day fasted versus 4-day fasted respectively. $\log_2$ fc values are also based on this formula. Outliers were excluded for this analysis (Fig \ref{fig:mahalanobis-outliers}).",
  '4vR': r"Differentially significant metabolites in 4-day fasted states which are similar in both cave populations. Cf. Table \ref{table:significant-30vR} with the difference that this table compares 4-day fasted versus refed conditions. Thus, metabolites displayed as upregulated in this table are differentially upregulated in 4-day fasted cave fish versus refed cave fish using surface fish as a baseline for comparison. Outliers are not included in this analysis.",
  '30v4': r"Differentially significant metabolites in 30-day vs 4-day fasting which are similar between both cave populations. Cf. Table \ref{table:significant-30vR} with the differences that this table compares 30-day fasted versus 4-day fasted conditions. Thus, metabolites displayed as upregulated in this table are differentially upregulated in 30-day fasted cave fish versus 4-day fasted cave fish using surface fish as a baseline for comparison. Outliers are not included in this analysis.",
}

labels = {
  '30vR': 'significant-30vR',
  '4vR':  'significant-4vR',
  '30v4':  'significant-30v4',
}

for comp in comparisons:
    #https://tex.stackexchange.com/questions/331716/newline-in-multirow-environment
    tabletex = r'''
\begin{table*}[t]
\centering
\caption{
{\bf $caption } }
\begin{tabular}
{p{1.5cm} | p{3cm} r|p{3cm} r|p{3cm} r}
\hline
{\bf Category} & {\bf Brain} & {$\log_2\mathrm{fc}$} & {\bf Muscle} & {$\log_2\mathrm{fc}$} & {\bf Liver} & {$\log_2\mathrm{fc}$} \\ \hline
$lines
\end{tabular}
\begin{flushleft}
$desc
\end{flushleft}
\label{table:$label}
\end{table*}
    '''.replace(
      '$lines','\n\\hline\n'.join(catlines[comp,category] for category in categories.values())).replace(
      '$desc',descs[comp]
      ).replace(
      '$caption',
      captions[comp]
      ).replace(
      '$label',
      labels[comp]
      )
    #print(f'\n\n\n% {comp}\n\n\n')
    if args.out_dir is not None:
        with open(os.path.join(args.out_dir,f'{comp}.tex'),'w') as f:
            f.write(tabletex)
    else:
        print(tabletex)
