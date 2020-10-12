from argparse import ArgumentParser
from pandas import read_csv
from cavefinomics import AstyanaxLi
from seaborn import heatmap, distplot
from seaborn import catplot, FacetGrid
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import matplotlib
from numpy import log10, isfinite, any, all, quantile, nan, isnan
from pandas import DataFrame, concat, melt, Categorical
import pathlib
import os
from cavefinomics import run_opls
from sklearn.preprocessing import StandardScaler, Normalizer
from numpy import log10, mean, isnan, quantile, nan, array
import matplotlib.patches as patches

parser = ArgumentParser(description="Lipid category table.")
parser.add_argument("--input-dir", type=str, help="Input directory.")
parser.add_argument("--outlier", type=str, help="Outlier directory name.")
parser.add_argument('--output-dir', type=str, help="Output directory.")
args = parser.parse_args()


#https://tex.stackexchange.com/questions/112343/beautiful-table-samples

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
polarities = ['positive','negative']
conditions = {'4d Starved':'4d', '30d Starved':'30d', 'Refed':'Ref'}
conditions_short = ['4d\nStarved', '30d\nStarved', 'Refed']
conditions_really_short = ['4', r'30', 'R']
comparisons = {'PvS':('Pachon','Surface'), 'TvS':('Tinaja','Surface'), 'PvT':('Pachon','Tinaja')}

cats = {
  'Fatty.Acyls':'Fatty Acyls',
  'Glycerolipids':'Glycerolipids',
  'Glycerophospholipids':'Glycerophospholipids',
  'Sphingolipids':'Sphingolipids',
  'Sterol.Lipids':'Sterol Lipids',
}
classes = {
  'Ceramides':'Ceramides',
  'Diradylglycerols':'Diradyl-glycerols',
  'Docosanoids':'Docosanoids',
  'Fatty.Acids.and.Conjugates':'Fatty Acids',
  'Fatty.esters':'Fatty esters',
  'Glycerophosphocholines':'Glycerophospho-cholines',
  'Glycerophosphoethanolamines':'Glycerophospho-ethanolamines',
  'Glycerophosphoglycerols':'Glycerophospho-glycerols',
  'Glycerophosphoglycerophosphoglycerols':'Glycerophospho-glycerophospho-glycerols',
  'Glycerophosphoinositols':'Glycerophospho-inositols',
  'Glycerophosphoserines':'Glycerophospho-serines',
  'Glycosyldiradylglycerols':'Glycosyldiradyl-glycerols',
  'Monoradylglycerols':'Monoradylglycerols',
  'Neutral.glycosphingolipids':'Neutral glycosphingolipids',
  'Oxidized.glycerophospholipids':'Oxidized glycerophospho-lipids',
  'Phosphosphingolipids':'Phosphosphingo-lipids',
  'Sphingoid.bases':'Sphingoid bases',
  'Sterols':'Sterols',
  'Triradylglycerols':'Triradylglycerols',
}

cattypes = {'Classes':classes,'Categories':cats}

sig_table = []
for tissue in tissues:
    for condition in conditions.values():
        for cattype in cattypes:
            for comparison in comparisons:
                d = read_csv(os.path.join(args.input_dir,f'{args.outlier}/{tissue}/{condition}/{cattype}/{comparison}.csv'))
                d['Tissue'] = tissue
                d['Condition'] = condition
                d['CatType'] = cattype
                d['Comparison'] = comparison
                sig_table.append(d)

sig_table = concat(sig_table)
cols = list(sig_table.columns)
cols[0] = 'Name'
sig_table.columns = cols
sig_table = sig_table.set_index('Name')
sig_table['Tissue'] = Categorical(sig_table['Tissue'], tissues)
sig_table['Condition'] = Categorical(sig_table['Condition'], conditions.values())
sig_table['Comparison'] = Categorical(sig_table['Comparison'], comparisons)


def make_comp_text(comp):
    return ' vs. '.join(comp)

def bold(u):
    return f'\textbf{{{str(u)}}}'

table_cols = '| m{3cm} | ' + ' || '.join([' | '.join(['m{0.225cm}']*9)]*3) + ' | '
table_comparison_header = r'\multicolumn{1}{c|}{} & ' + ' & '.join([f'\\multicolumn{{9}}{{c|}}{{{bold(make_comp_text(comp))}}}' for comp in comparisons.values()]) + r' \\'
table_tissue_header = r'\multicolumn{1}{c|}{} & ' + ' & '.join([f'\\multicolumn{{3}}{{c|}}{{{bold(tissue)}}}' for tissue in tissues]*3) + r' \\'
table_condition_header = r'\multicolumn{1}{c|}{} & ' + ' & '.join([f'{bold(condition)}' for condition in conditions_really_short]*9) + r' \\'

def make_subtable_for_cat(cat):
    d = sig_table.loc[cat,['Comparison','Tissue','Condition','Pr(>|z|)','Estimate']]
    if len(d.index) < 27:
        # fill in missing values
        for comparison in comparisons:
            for tissue in tissues:
                for condition in conditions.values():
                    if not ((d['Comparison'] == comparison) & (d['Tissue'] == tissue) & (d['Condition'] == condition)).any():
                        new_row = DataFrame([{
                            'Name':cat,
                            'Comparison':comparison,
                            'Tissue':tissue,
                            'Condition':condition,
                            'Pr(>|z|)':nan,
                            'Estimate':nan}])
                        new_row = new_row.set_index('Name')
                        d = concat((d,new_row))
    return d.sort_values(['Comparison','Tissue','Condition'])

def get_sigs(level,cat):
    d = make_subtable_for_cat(cat)

    def conserved_highlight(input,comp,tissue,condition):
        if comp != 'PvS' and comp != 'TvS':
            return input
        u = float(d.loc[(d['Comparison'] == 'PvS') & (d['Tissue'] == tissue) & (d['Condition'] == condition),'Pr(>|z|)'])
        uup = float(d.loc[(d['Comparison'] == 'PvS') & (d['Tissue'] == tissue) & (d['Condition'] == condition),'Estimate']) > 0.
        v = float(d.loc[(d['Comparison'] == 'TvS') & (d['Tissue'] == tissue) & (d['Condition'] == condition),'Pr(>|z|)'])
        vup = float(d.loc[(d['Comparison'] == 'PvS') & (d['Tissue'] == tissue) & (d['Condition'] == condition),'Estimate']) > 0.
        if u < level and v < level and uup == vup:
            return '\cellcolor{cellhl}{'+input+'}'
        else:
            return input

    ps = d.loc[cat,'Pr(>|z|)']
    up = d.loc[cat,'Estimate'] > 0.
    cs = d.loc[cat,'Comparison']
    ts = d.loc[cat,'Tissue']
    cnds = d.loc[cat,'Condition']
    row = [conserved_highlight(r'$+$' if u else r'$-$', comp,tissue,condition) if p < level else (' ' if isfinite(p) else r' $\varnothing$ ') for p,u,comp,tissue,condition in zip(ps,up,cs,ts,cnds)]
    if len(row) == 27:
        return row
    else:
        return [' x ']*27


#https://tex.stackexchange.com/questions/130818/how-to-draw-a-double-hline-in-a-table-without-interrupting-vertical-lines
#https://tex.stackexchange.com/questions/152101/how-to-draw-two-hline-and-two-vertical-line

labels = {
  'Classes': 'table:sig-lipid-classes',
  'Categories': 'table:sig-lipid-categories',
}

captions = {
  'Classes': 'Interpopulation Differences in Abundance of Lipid Classes',
  'Categories': 'Interpopulation Differences in Abundance of Lipid Categories',
}

descriptions = {
  'Classes': r"Peak intensities for all lipids in a given class (determined from the LipidMaps ``\texttt{MAIN\_CLASS}'' attribute) were summed to yield a total intensity for each class which is either significantly (at the $p<0.05$ level) up-- ($+$) or down--regulated ($-$) in a given cave population with respect to surface (first two comparisons, top row) or the Pach\'{o}n cave population with respect to the Tinaja cave population (last comparison, top row). The sample set for each tissue / feeding state combination consists of six individuals from each population as shown in Fig \ref{fig:exp-setup}. P--values were obtained from the OPLS / GLM approach described in Methods. $\varnothing$ denotes classes which were not detected in a given sample set. LipidMaps also possesses a ``\texttt{CATEGORY}'' attribute that provides a more coarse--grained classification of lipid species, which is used as a basis for a similar analysis shown in Table \ref{table:sig-lipid-categories}.",
  'Categories': r"Significantly up-- or down--regulated lipid categories based on summed peak intensities for every lipid species belonging to a given LipidMaps ``\texttt{CATEGORY}'' label. Compare Table \ref{table:sig-lipid-classes} based on the ``\texttt{MAIN\_CLASS}'' attribute.",
}

level = 0.05
for cattype in cattypes:
    table_data = []
    for cat,name in cattypes[cattype].items():
        table_data.append(f'{name} & ' + ' & '.join(get_sigs(level,cat)) + r' \\')
        table_data.append(r'\hline')

    table_data = '\n'.join(table_data)

    tabletex = r'''
\begin{table*}[t]
\begin{adjustwidth}{-1cm}{}
\centering
\caption{
{\bf $caption } }
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
\begin{flushleft}
$desc
\end{flushleft}
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
          descriptions[cattype]
          ).replace(
          '$caption',
          captions[cattype]
          ).replace(
          '$label',
          labels[cattype]
          )

    print(tabletex)
    with open(os.path.join(args.output_dir,f'{cattype}.tex'),'w') as f:
        f.write(tabletex)
