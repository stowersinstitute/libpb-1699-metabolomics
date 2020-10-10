from argparse import ArgumentParser
import os

from numpy import array, log10, log2
from pandas import DataFrame, Series, concat
import matplotlib.pyplot as plt
import matplotlib
from seaborn import heatmap
import seaborn as sns
from functools import reduce
import math
from pprint import pprint
from collections import defaultdict
from pandas import read_csv

from cavefinomics import AstyanaxMe

#https://stackoverflow.com/questions/34706845/change-xticklabels-fontsize-of-seaborn-heatmap
sns.set(font_scale=1.4)

parser = ArgumentParser(description="Heatmap of conserved metabolites.")
parser.add_argument("--output-dir", type=str, help="Output directory.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
args = parser.parse_args()

conditions = {'4d':'4d Starved', '30d':'30d Starved', 'Ref':'Refed'}
pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
comparisons = {'30vR':('30d Starved','Refed'),'4vR':('4d Starved','Refed'),'30v4':('30d Starved','4d Starved')}
categories = {"Aminoacids":'Amino acids',"Carbohydrates_-CCM": 'Carbohydrates / CCM',"Fattyacids":'Fatty acids',"Misc._-_sec.metabolites":'Misc',"Nucleotides":'Nucleotides'}
if args.exclude_outlier:
    outlier_text = 'kein-Ausreißern'
else:
    outlier_text = 'mit-Ausreißern'

metabolites = {}
for tissue in tissues:
    for comp,groups in comparisons.items():
        sig = []
        for cat,category in categories.items():
            data = read_csv(f"out/work/primary/glm/singlefactor/{outlier_text}/{cat}/{tissue}/CvS/{comp}.csv").rename({'Pr(>|z|)':'p'},axis=1)
            cols = list(data.columns)
            cols[0] = 'Name'
            data.columns = cols
            data = data.set_index('Name')
            data['Tissue'] = tissue
            data['Comparison'] = comp
            sig.append(data)
        sig = concat(sig)
        metabolites[tissue,comp] = sig



remap_cond = {
  '4d Starved': '4d Fasted',
  '30d Starved': '30d Fasted',
  }

for comp,groups in comparisons.items():
    gridspec_kw = {"height_ratios":[1.], "width_ratios" : [3.,3.,3.]}
    fig,ax = plt.subplots(nrows=1,ncols=3,figsize=(12, 8), gridspec_kw=gridspec_kw)
    fig.suptitle(' vs '.join((remap_cond[g] if g in remap_cond else g for g in groups)),fontsize='xx-large',fontweight='bold')

    vmin = -1.7
    vmax = 1.7

    for j,tissue in zip(range(3),tissues):
        d = metabolites[tissue,comp].sort_values('p')
        reg = -d['p'].apply(log10)
        reg.loc[d['Estimate'] < 0] = -reg
        reg.index = [u if v > 0.05 else ''.join((u,'*')) for u,v in zip(reg.index,d['p'])]
        N = 10
        reg = reg.iloc[:N].sort_values(ascending=False)

        heatmap(array([list(reg.iloc[:N].values)]).T, vmin=vmin, vmax=vmax, annot=array([list(reg.iloc[:N].index)]).T, fmt = '', ax=ax[j], cbar=False, xticklabels=False, yticklabels=j==0, cmap='coolwarm')
        ax[j].set_yticks([], minor=[])
        ax[j].set_xticks([], minor=[])
        ax[j].set_title(tissue,fontsize='x-large',fontweight='bold')

    filename=f'{args.output_dir}/{outlier_text}/primary-shared-starvation-response-{comp}.pdf'
    plt.savefig(filename,bbox_inches='tight',transparent=True,pad_inches=0)

fig,ax = plt.subplots(figsize=(3, 12))
# draw c bar
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
sm = matplotlib.cm.ScalarMappable(cmap='coolwarm', norm=norm)
sm.set_array([])
cax = fig.add_axes([0.1,0.1,0.25,0.8])
fig.colorbar(sm, cax=cax)
#cax.set_title(r'$-d\log_{10}p$',fontsize='xx-large')
cax.yaxis.set_major_locator(plt.FixedLocator([vmin, 0., vmax]))
cax.set_yticklabels([f'$p = {10.**vmin:.2f}$',r'$p = 1$',f'$p = {10.**(-vmax):.2f}$'],fontsize='xx-large',fontweight='bold')
cax.text(0.5,1.1,'Up in \ncave+fasted',size=20,ha='center',va='center',fontweight='bold',transform=ax.transAxes)
cax.text(0.5,-0.1,'Down in \ncave+fasted',size=20,ha='center',va='center',fontweight='bold',transform=ax.transAxes)
# clear extra axes
# hide graphics
ax.set_yticks([], minor=[])
ax.set_xticks([], minor=[])
ax.patch.set_visible(False)
for s in ["top", "bottom", "left", "right"]:
    ax.spines[s].set_visible(False)
filename=f'{args.output_dir}/{outlier_text}/primary-shared-starvation-response-colorbar.pdf'
plt.savefig(filename)
