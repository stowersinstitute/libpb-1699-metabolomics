from argparse import ArgumentParser
import os

from numpy import array, log10, log2
from pandas import DataFrame, Series, concat
import matplotlib.pyplot as plt
import matplotlib
from seaborn import heatmap
from functools import reduce
import math
from pprint import pprint
from collections import defaultdict
from pandas import read_csv

from cavefinomics import AstyanaxMe, MaMammalianMe

conditions = {'4d':'4d Starved', '30d':'30d Starved', 'Ref':'Refed'}
pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
comparisons = {'30vR':('30d Starved','Refed'),'4vR':('4d Starved','Refed'),'30v4':('30d Starved','4d Starved')}
#comparisons = {'30vR':('30d Starved','Refed')}
categories = {"Aminoacids":'Amino acids',"Carbohydrates_-CCM": 'Carbohydrates / CCM',"Fattyacids":'Fatty acids',"Misc._-_sec.metabolites":'Misc',"Nucleotides":'Nucleotides'}
outlier = 'no-outliers'

metabolites = {}
for tissue in tissues:
    for comp,groups in comparisons.items():
        sig = []
        for cat,category in categories.items():
            data = read_csv(f"out/work/primary/glm/singlefactor/{outlier}/{cat}/{tissue}/CvS/{comp}.csv").rename({'Pr(>|z|)':'p'},axis=1)
            cols = list(data.columns)
            cols[0] = 'Name'
            data.columns = cols
            #if comp == 'PvT':
                #data['Estimate'] = -data['Estimate']
            data = data.set_index('Name')
            data['Tissue'] = tissue
            data['Comparison'] = comp
            sig.append(data)
            #print(data)
        sig = concat(sig)
        #print(sig)
        metabolites[tissue,comp] = sig

#print(metabolites)

for comp,groups in comparisons.items():
    gridspec_kw = {"height_ratios":[1.], "width_ratios" : [3.,3.,3.,1.,1]}
    fig,ax = plt.subplots(nrows=1,ncols=5,figsize=(12, 8), gridspec_kw=gridspec_kw)
    fig.suptitle(' vs '.join(groups))

    vmin = -1.7
    vmax = 1.7

    for j,tissue in zip(range(3),tissues):

        d = metabolites[tissue,comp].sort_values('p')
        #print(d)
        reg = -d['p'].apply(log10)
        reg.loc[d['Estimate'] < 0] = -reg
        reg.index = [u if v > 0.05 else ''.join((u,'*')) for u,v in zip(reg.index,d['p'])]
        reg = reg.iloc[:20].sort_values(ascending=False)
        #print(reg)

        #print(reg.iloc[:20])
        #print(list(reg.iloc[:20].values))
        #print(list(reg.iloc[:20].index))
        #print(array([list(reg.iloc[:20].index)]).shape)
        heatmap(array([list(reg.iloc[:20].values)]).T, vmin=vmin, vmax=vmax, annot=array([list(reg.iloc[:20].index)]).T, fmt = '', ax=ax[j], cbar=False, xticklabels=False, yticklabels=j==0, cmap='coolwarm')
        ax[j].set_yticks([], minor=[])
        ax[j].set_xticks([], minor=[])
        ax[j].set_title(tissue)

    # clear extra axes
    # hide graphics
    ax[3].set_yticks([], minor=[])
    ax[3].set_xticks([], minor=[])
    ax[3].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        ax[3].spines[s].set_visible(False)
    # hide graphics
    ax[4].set_yticks([], minor=[])
    ax[4].set_xticks([], minor=[])
    ax[4].patch.set_visible(False)
    for s in ["top", "bottom", "left", "right"]:
        ax[4].spines[s].set_visible(False)

    # draw c bar
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = matplotlib.cm.ScalarMappable(cmap='coolwarm', norm=norm)
    sm.set_array([])
    cax = fig.add_axes([0.9,0.1,0.03,0.8])
    fig.colorbar(sm, cax=cax)
    cax.set_title(r'$-d\log_{10}p$')
    cax.yaxis.set_major_locator(plt.FixedLocator([vmin, vmax]))
    cax.set_yticklabels([f'{vmin:.1f}',f'{vmax:.1f}'])
    filename=f'/tmp/primary-shared-heatmap-{comp}.pdf'
    plt.savefig(filename)
