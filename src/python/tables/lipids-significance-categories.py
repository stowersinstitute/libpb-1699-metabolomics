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
#http://cpansearch.perl.org/src/LIMAONE/LaTeX-Table-v1.0.6/examples/examples.pdf
#https://tex.stackexchange.com/questions/113101/how-can-i-avoid-the-gap-between-cellcolor-colored-cells-created-by-tabucline

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
    return f'\\textbf{{{str(u)}}}'

table_cols = '| r{3cm} | ' + ' | '.join([' | '.join(['l']*9)]*3) + ' | '
table_comparison_header = r'\multicolumn{1}{c}{} & ' + ' & '.join([f'\\multicolumn{{9}}{{c}}{{{bold(make_comp_text(comp))}}}' for comp in comparisons.values()]) + r' \\'
table_tissue_header = r'\multicolumn{1}{c}{} & ' + ' & '.join([f'\\multicolumn{{3}}{{c}}{{{bold(tissue)}}}' for tissue in tissues]*3) + r' \\'
table_condition_header = r'\multicolumn{1}{c}{} & ' + ' & '.join([f'\\multicolumn{{1}}{{c}}{{{bold(condition)}}}' for condition in conditions_really_short]*9) + r' \\'

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
    d['Tissue'] = Categorical(d['Tissue'], tissues)
    d['Condition'] = Categorical(d['Condition'], conditions.values())
    d['Comparison'] = Categorical(d['Comparison'], comparisons)
    return d.sort_values(['Comparison','Tissue','Condition'])

def get_sigs(level,cat):
    d = make_subtable_for_cat(cat)

    def conserved_highlight(input,comp,tissue,condition,ismuscle):
        if comp != 'PvS' and comp != 'TvS':
            if ismuscle:
                return '\cellcolor{structcell}{'+input+'}'
            else:
                return input
        u = float(d.loc[(d['Comparison'] == 'PvS') & (d['Tissue'] == tissue) & (d['Condition'] == condition),'Pr(>|z|)'])
        uup = float(d.loc[(d['Comparison'] == 'PvS') & (d['Tissue'] == tissue) & (d['Condition'] == condition),'Estimate'])
        v = float(d.loc[(d['Comparison'] == 'TvS') & (d['Tissue'] == tissue) & (d['Condition'] == condition),'Pr(>|z|)'])
        vup = float(d.loc[(d['Comparison'] == 'TvS') & (d['Tissue'] == tissue) & (d['Condition'] == condition),'Estimate'])
        if u < level and v < level and uup*vup > 0.:
            return '\cellcolor{cellhl}{'+input+'}'
        else:
            if ismuscle:
                return '\cellcolor{structcell}{'+input+'}'
            else:
                return input

    ps = list(d.loc[cat,'Pr(>|z|)'])
    up = d.loc[cat,'Estimate'] > 0.
    ismuscle = list(d.loc[cat,'Tissue'] == 'Muscle')
    cs = d.loc[cat,'Comparison']
    ts = d.loc[cat,'Tissue']
    cnds = d.loc[cat,'Condition']
    row = [
      conserved_highlight(
          (r'$+$' if u else r'$-$') if p < level else (' ' if isfinite(p) else r' $\varnothing$ '),
          comp,tissue,condition,mus) for p,u,comp,tissue,condition,mus in zip(ps,up,cs,ts,cnds,ismuscle)]
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

# square
#https://tex.stackexchange.com/questions/369429/write-block-character-with-char/369568#369568

descriptions = {
  'Classes': r"Significant increased ($+$) or decreased ($-$) lipid classes based on summed peak intensities for every lipid species belonging to a given LipidMaps ``\texttt{MAIN\_CLASS}'' label. Compare Table \ref{table:sig-lipid-classes} based on the ``\texttt{MAIN\_CLASS}'' attribute. $\varnothing$ denotes classes which were not detected in a given sample set.",

  'Categories': r"Peak intensities for all lipids in a given category (determined from the LipidMaps ``\texttt{CATEGORY}'' attribute) were summed to yield a total intensity for each category which is either significantly (at the $p<0.05$ level) up-- ($+$) or down--regulated ($-$) in a given cave population with respect to surface (Pach\'{o}n versus surface and Tinaja versus surface, top row) or the Pach\'{o}n cave population with respect to the Tinaja cave population (last comparison, top row). The sample set for each tissue / feeding state combination consists of six individuals from each population as shown in Fig \ref{fig:exp-setup}. P--values were obtained from the OPLS / GLM approach described in Methods. Coloring (\textcolor{cellhl}{$\blacksquare\!$}) indicates a class that agrees in significance and directionality between both cave populations and is thus may be related to cave adaptation. LipidMaps also possesses a ``\texttt{CATEGORY}'' attribute that provides a more coarse--grained classification of lipid species, which is used as a basis for a similar analysis shown in Table \ref{table:sig-lipid-categories}.",
}

level = 0.05
for cattype in cattypes:
    table_data = []
    for cat,name in cattypes[cattype].items():
        table_data.append(f'\multicolumn{{1}}{{r|}}{{{name}}} & ' + ' & '.join(get_sigs(level,cat)) + r' \\')
        table_data.append(r'\cline{2-28}')

    table_data = '\n'.join(table_data)

    tabletex = r'''
\begin{table*}[t]
\begin{adjustwidth}{-1cm}{}
\centering
\resizebox{1.05\textwidth}{!}{
\begin{tabular}
{$cols}
$table_comparison_header
$table_tissue_header
$table_condition_header
\cline{2-28}
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
          descriptions[cattype]
          ).replace(
          '$caption',
          captions[cattype]
          ).replace(
          '$label',
          labels[cattype]
          )

    with open(os.path.join(args.output_dir,f'{cattype}.tex'),'w') as f:
        f.write(tabletex)
