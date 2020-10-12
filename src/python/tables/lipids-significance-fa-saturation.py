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
#plt.style.use('seaborn')
#plt.style.use('Solarize_Light2')
#plt.style.use('fivethirtyeight')
#matplotlib.rc('image', cmap='Pastel2')
#plt.set_cmap('jet')
#plt.set_cmap('Pastel2')
#https://colorhunt.co/palette/192017
#matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#99b898", "#feceab", "#ff847c","#e84a5f"])
#https://colorhunt.co/palette/191947
#matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#111d5e", "#c70039", "#f37121","#ffbd69"])
#https://colorhunt.co/palette/192164
#matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#e7305b", "#e2979c", "#f7f5dd","#9bdeac"])
#https://colorhunt.co/palette/189676
#matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["#726a95", "#709fb0", "#a0c1b8","#f4ebc1"])
#https://colorhunt.co/palette/202236
# beige
#https://colorhunt.co/palette/202295
# pink
#https://colorhunt.co/palette/195720
# blue/orange
#https://colorhunt.co/palette/195946
# blue/pink
#https://colorhunt.co/palette/196113
# blues
#https://colorhunt.co/palette/195717
# desert camo

parser = ArgumentParser(description="Lipid category table.")
parser.add_argument("--input-dir", type=str, help="Input directory.")
parser.add_argument("--outlier", type=str, help="Outlier directory name.")
parser.add_argument('--output-dir', type=str, help="Output directory.")
args = parser.parse_args()

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
  'Monounsaturated.Fatty.Acids': 'Monounsaturated Fatty Acids',
  'Polyunsaturated.Fatty.Acids': 'Polyunsaturated Fatty Acids',
  'Saturated.Fatty.Acids': 'Saturated Fatty Acids',
}
wanted_classes = {
  'Monounsaturated.Fatty.Acids': 'Monounsat. FAs',
  'Polyunsaturated.Fatty.Acids': 'Polyunsat. FAs',
  'Saturated.Fatty.Acids': 'Saturated FAs',
}

cattypes = {'Classes':classes,'Categories':cats}

sig_table = []
for tissue in tissues:
    for condition in conditions.values():
        for comparison in comparisons:
            d = read_csv(os.path.join(args.input_dir,f'{args.outlier}/{tissue}/{condition}/fas/{comparison}.csv'))
            d['Tissue'] = tissue
            d['Condition'] = condition
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
#print(sig_table)


def make_comp_text(comp):
    return ' vs. '.join(comp)

table_cols = '| r{3cm} | ' + ' | '.join([' | '.join(['l']*9)]*3) + ' | '
table_comparison_header = r'\multicolumn{1}{c}{} & ' + ' & '.join([f'\\multicolumn{{9}}{{c|}}{{{make_comp_text(comp)}}}' for comp in comparisons.values()]) + r' \\'
table_tissue_header = r'\multicolumn{1}{c}{} & ' + ' & '.join([f'\\multicolumn{{3}}{{c|}}{{{tissue}}}' for tissue in tissues]*3) + r' \\'
table_condition_header = r'\multicolumn{1}{c}{} & ' + ' & '.join([f'{condition}' for condition in conditions_really_short]*9) + r' \\'

#https://tex.stackexchange.com/questions/178622/how-to-place-a-line-for-n-a-at-a-cell-of-a-table
# row height
#https://tex.stackexchange.com/questions/159257/increase-latex-table-row-height

def get_sigs(level,cat):
    d = sig_table.loc[cat,['Comparison','Tissue','Condition','Pr(>|z|)','Estimate']]
    if len(d.index) < 27:
        for comparison in comparisons:
            for tissue in tissues:
                for condition in conditions.values():
                    if not ((d['Comparison'] == comparison) & (d['Tissue'] == tissue) & (d['Condition'] == condition)).any():
                        pass
                        new_row = DataFrame([{'Name':cat,'Comparison':comparison,'Tissue':tissue,'Condition':condition,'Pr(>|z|)':nan,'Estimate':nan}])
                        new_row = new_row.set_index('Name')
                        d = concat((d,new_row))
    d = d.sort_values(['Comparison','Tissue','Condition'])
    ps = d.loc[cat,'Pr(>|z|)']
    up = d.loc[cat,'Estimate'] > 0.
    ismuscle = list(d.loc[cat,'Tissue'] == 'Muscle')
    conserved = [True if ps.iloc[k] < level and ps.iloc[k+9] < level and up.iloc[k] == up.iloc[k+9] else False for k in range(9)]
    conserved = conserved+conserved+[False]*9
    def make_symbol(p,u,conserved):
        def add_cell_color(x):
            if conserved:
                return r'\cellcolor{cellhl}{'+x+'}'
            else:
                return x
        if p < level:
            if u:
                return add_cell_color(r'$+$')
            else:#else
                return add_cell_color(r'$-$')
        elif isfinite(p):
            return ' '
        else:
            return ' --- '
    def highlight_muscle(input, ismuscle):
        if ismuscle:
            return '\cellcolor{structcell}{'+input+'}'
        else:
            return input
    row = [highlight_muscle(make_symbol(p, u, c ),mus) for p,u,mus,c in zip(ps,up,ismuscle,conserved)]
    if len(row) == 27:
        return row
    else:
        return [' x ']*27


#https://tex.stackexchange.com/questions/130818/how-to-draw-a-double-hline-in-a-table-without-interrupting-vertical-lines
#https://tex.stackexchange.com/questions/152101/how-to-draw-two-hline-and-two-vertical-line

label = 'table:sig-lipid-categories-fas'

caption = 'Intra--population Differences in Fatty Acid Saturation'

description = r"For lipids corresponding to free fatty acids, sauturation was calculated based on the presence of double bonds in LipidMaps structural data and used to classify each LMID as either saturated, monounsaturated, or polyunsaturated. Significance values were again calculated using an OPLS / Bayesian GLM workflow. Coloring and markings as before."

level = 0.05

table_data = []
for cat,name in wanted_classes.items():
    table_data.append(f'\\multicolumn{{1}}{{r|}}{{{name}}} & ' + ' & '.join(get_sigs(level,cat)) + r' \\')
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
      description
      ).replace(
      '$caption',
      caption
      ).replace(
      '$label',
      label
      )

with open(os.path.join(args.output_dir,f'fas.tex'),'w') as f:
    f.write(tabletex)
