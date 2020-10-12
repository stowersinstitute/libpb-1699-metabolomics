from argparse import ArgumentParser
import matplotlib.pyplot as plt
from numpy import histogram,max,argmax,nan,isnan
from pandas import read_csv
from pandas import concat
from cavefinomics import AstyanaxLi
import json
import os

parser = ArgumentParser(description="Compute Inchi keys for compounds.")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument(
    "--lipidmaps-json", type=str, help="LipidMaps JOSN file.",
)
parser.add_argument("--input-dir", type=str, help="Input directory.")
parser.add_argument("--outlier", type=str, help="Outlier directory name.")
parser.add_argument('--output', type=str, help="Output.")
args = parser.parse_args()

ali = AstyanaxLi(
  lipids_normalized=args.lipids_normalized,
  lipidmaps_js=args.lipidmaps_json,
  )

with open(args.lipidmaps_json) as f:
    lm = json.load(f)

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = ['4d Starved', '30d Starved', 'Refed']
conditions_really_short = ['4', r'30', 'R']
condmap = {'30d':'30d Starved', '4d':'4d Starved', 'Ref': 'Refed'}
comparisons = {'PvS':('Pachon','Surface'), 'TvS':('Tinaja','Surface'), 'PvT':('Pachon','Tinaja')}


def get_lipidmap_id(inchikey):
    if inchikey in ali.lipidmaps_inchikey:
        return ali.lipidmaps_inchikey[inchikey]
    elif '-'.join(inchikey.split('-')[:2]) in ali.lipidmaps_inchikey2:
        return ali.lipidmaps_inchikey2['-'.join(inchikey.split('-')[:2])]
    elif inchikey.split('-')[0] in ali.lipidmaps_inchikey1:
        return ali.lipidmaps_inchikey1[inchikey.split('-')[0]]
    return nan

def get_lipidmap_name(lmid):
    if 'NAME' in lm[lmid]:
        return lm[lmid]['NAME']
    else:
        return nan

def is_name_common(name):
    if len(name) >= 3 and name[2] == '(':
        return False
    elif len(name) >= 4 and name.startswith('Cer('):
        return False
    elif len(name) >= 7 and name[3:].startswith('Cer('):
        return False
    elif name == '15:1(9Z)':
        return False
    elif name == 'MGDG(16:0/18:1(9Z))':
        return False
    else:
        return True

wanted_lmids = {
}

def rename(name):
    namemap = {
      '9Z,12Z,15Z,18Z,21Z-tetracosapentaenoic acid': 'tetracosapentaenoic acid',
      #'alpha-Linolenic acid': 'Linolenic acid'
      }
    if name in namemap:
        return namemap[name]
    else:
        return name

datasets = []
for tissue in tissues:
    for cond in condmap:
        for comp in comparisons:
            d = read_csv(os.path.join(args.input_dir,f'{args.outlier}/{tissue}/{cond}/{comp}.csv'),index_col=0)
            d['Tissue'] = tissue
            d['Condition'] = condmap[cond]
            d['Comparison'] = comp
            d['LMID'] = [get_lipidmap_id(k) for k in d.index]
            d['Name'] = [get_lipidmap_name(l) for l in d['LMID']]
            datasets.append(d)
data = concat(datasets,axis=0).dropna()
n_lmids = len(set(d['LMID']))
#print(f'Num unique LMIDs: {n_lmids}')
data = data.loc[[is_name_common(name) or lmid in wanted_lmids for name,lmid in zip(data['Name'],data['LMID'])],:]
data['Name'] = [rename(name) for name in data['Name']]
#print(f'Num with human-readable names: {len(data.index)}')

table_data = []
sig_data = {}
for name in sorted(set(data['Name'])):
    for tissue in tissues:
        for condition in condmap.values():
            for comp in comparisons:
                est = float(data.loc[(data['Tissue'] == tissue) & (data['Condition'] == condition) & (data['Comparison'] == comp) & (data['Name'] == name),'Estimate'])
                p = float(data.loc[(data['Tissue'] == tissue) & (data['Condition'] == condition) & (data['Comparison'] == comp) & (data['Name'] == name),'Pr(>|z|)'])
                if p < 0.05:
                    if est > 0:
                        sig_data[tissue,condition,comp,name] = 1
                    else:
                        sig_data[tissue,condition,comp,name] = -1
                else:
                    sig_data[tissue,condition,comp,name] = 0
def conserved(tissue,condition,comp,name):
    if comp == 'PvS' or comp == 'TvS':
        return sig_data[tissue,condition,'PvS',name] == sig_data[tissue,condition,'TvS',name]
    else:
        return False

def highlight_muscle(input, ismuscle):
    if ismuscle:
        return '\cellcolor{structcell}{'+input+'}'
    else:
        return input

def get_sigs(name):
    r = []
    for comp in comparisons:
        for tissue in tissues:
            for condition in condmap.values():
                s = sig_data[tissue,condition,comp,name]
                if s != 0:
                    if s > 0:
                        if conserved(tissue,condition,comp,name):
                            r.append(r'\cellcolor{cellhl}{$+$}')
                        else:
                            r.append(highlight_muscle(r'$+$', tissue=='Muscle'))
                    else:
                        if conserved(tissue,condition,comp,name):
                            r.append(r'\cellcolor{cellhl}{$-$}')
                        else:
                            r.append(highlight_muscle(r'$-$',tissue=='Muscle'))
                else:
                    r.append(highlight_muscle('', tissue=='Muscle'))
    return r
for name in sorted(set(data['Name'])):
    table_data.append(f'{name} & ' + ' & '.join(get_sigs(name)) + r' \\\hline')
table_data = '\n'.join(table_data)

def make_comp_text(comp):
    return ' vs. '.join(comp)

table_cols = '| m{3cm} | ' + ' || '.join([' | '.join(['m{0.225cm}']*9)]*3) + ' | '
table_comparison_header = r'\multicolumn{1}{c|}{} & ' + ' & '.join([f'\\multicolumn{{9}}{{c|}}{{{make_comp_text(comp)}}}' for comp in comparisons.values()]) + r' \\'
table_tissue_header = r'\multicolumn{1}{c|}{} & ' + ' & '.join([f'\\multicolumn{{3}}{{c|}}{{{tissue}}}' for tissue in tissues]*3) + r' \\'
table_condition_header = r'\multicolumn{1}{c|}{} & ' + ' & '.join([f'{condition}' for condition in conditions_really_short]*9) + r' \\'

tabletex = r'''
\begin{table*}[t]
\begin{adjustwidth}{-1cm}{}
\centering
\resizebox{1.05\textwidth}{!}{
\begin{tabular}
{$cols}
\cline{2-28}
$table_comparison_header
\cline{2-28}
$table_tissue_header
\cline{2-28}
$table_condition_header
\hline
$data
\end{tabular}
}
\vspace{0.75em}
\caption{
{\bf $caption }
$desc
}
\label{$label}
\end{adjustwidth}
\end{table*}
    '''.replace(
      '$cols',
      table_cols
      ).replace(
      '$table_comparison_header',
      table_comparison_header
      ).replace(
      '$table_tissue_header',
      table_tissue_header
      ).replace(
      '$table_condition_header',
      table_condition_header
      ).replace(
      '$data',
      table_data
      ).replace(
      '$desc',
      r'Lipids with common names were selected from the set of 447 identified lipids and analyzed using the O--PLS / GLM scheme. Many lipids exhibit a strongly conserved pattern between Pach\'{o}n and Tinaja, and these are highlighted in blue. Of note, omega--3 fatty acids appear to exhibit lower abundance in the 4-day fasted state in liver. Some metabolites, such as palmitate and oleic acid, overlap with primary metabolomics data.'
      ).replace(
      '$caption',
      'Common lipids'
      ).replace(
      '$label',
      'table:human-readable-lipids'
      )

with open(args.output,'w') as f:
    f.write(tabletex)
