from argparse import ArgumentParser
#from cavefinomics import AstyanaxMe
from cavefinomics import AstyanaxLi
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import matplotlib.ticker as plticker
from itertools import repeat, chain
from collections import defaultdict
flatten = chain.from_iterable
import os
import json
from pandas import concat, DataFrame, read_csv
from numpy import nan

parser = ArgumentParser(description="Lipids merged.")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument("--lipidmaps-json", type=str, help="Lipidmaps JSON.")
parser.add_argument("--lipidmaps-fa", type=str, help="Lipidmaps fa.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--out-cross-pop", type=str, help="Output (cross-pop)")
parser.add_argument("--out-cross-cond", type=str, help="Output (cross-cond)")
args = parser.parse_args()

with open(args.lipidmaps_fa) as f:
  lmfa = json.load(f)

with open(args.lipidmaps_json) as f:
    lm = json.load(f)

ali = AstyanaxLi(
    lipids_normalized=args.lipids_normalized,
    lipidmaps_js=args.lipidmaps_json,
    )

if args.exclude_outlier:
    outlier = 'kein-Ausreißern'
else:
    outlier = 'mit-Ausreißern'

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = ['4d Starved', '30d Starved', 'Refed']
comparisons = {'PvS':('Pachon','Surface'), 'TvS':('Tinaja','Surface'), 'PvT':('Pachon','Tinaja')}
condmap = {'30d':'30d Starved', '4d':'4d Starved', 'Ref': 'Refed'}

def read_sig_dataset(filepath,tissue,cond,comp):
    d = read_csv(filepath,index_col=0)
    d['Tissue'] = tissue
    d['Condition'] = cond
    d['Comparison'] = comp
    return d

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

def assignfa(lmid):
    if isinstance(lmid,str):
        mc = ali.lipidmaps[lmid]['MAIN_CLASS']
        if mc != "Fatty Acids and Conjugates [FA01]":
            return nan
        else:
            fa = lmfa[lmid]
            if 'unsat' in fa:
                if fa['unsat'] > 1:
                    return 'Polyunsaturated Fatty Acids'
                elif fa['unsat'] == 1:
                    return 'Monounsaturated Fatty Acids'
            elif 'sat' in fa:
                return 'Saturated Fatty Acids'
            return nan

cross_pop_significance = []
# cross population comparison
for tissue in tissues:
    for cond in condmap:
        for comp in comparisons:
            glm = read_sig_dataset(f'out/work/lipids/glm/singlefactor/{outlier}/{tissue}/{cond}/{comp}.csv',tissue,cond,comp)
            opls = read_sig_dataset(f'out/work/lipids/opls/{outlier}/{tissue}/{cond}/{comp}.csv',tissue,cond,comp)
            zscore = read_sig_dataset(f'out/work/lipids/zscore/{outlier}/{tissue}/{cond}/{comp}.csv',tissue,cond,comp)
            glm['LMID'] = glm.apply(lambda u: get_lipidmap_id(u.name),axis=1)
            glm['Category'] = glm.apply(lambda u: ali.lipidmaps[u['LMID']]['CATEGORY'] if isinstance(u['LMID'],str) else nan,axis=1)
            glm['MainClass'] = glm.apply(lambda u: ali.lipidmaps[u['LMID']]['MAIN_CLASS'] if isinstance(u['LMID'],str) else nan,axis=1)
            glm['Saturation'] = glm.apply(lambda u: assignfa(u['LMID']) if isinstance(u['LMID'],str) else nan,axis=1)
            cross_pop_significance.append(glm)
cross_pop_significance = concat(cross_pop_significance,axis=0).dropna(subset=['LMID'])
cross_pop_significance = cross_pop_significance.rename({'Pr(>|z|)':'p-val','Estimate':'Slope'},axis=1)
cross_pop_significance.index.name = 'Name'
if args.out_cross_pop:
    cross_pop_significance.to_csv(args.out_cross_pop)

def read_sig_dataset2(filepath,tissue,pop,comp):
    d = read_csv(filepath,index_col=0)
    d['Tissue'] = tissue
    d['Population'] = pop
    d['Comparison'] = comp
    return d
comparisons2 = ['30vR','4vR','30v4']

cross_cond_significance = []
# cross condition comparison
for tissue in tissues:
    for pop in pops:
        for comp in comparisons2:
            glm = read_sig_dataset2(f'out/work/lipids/glm/singlefactor/{outlier}/{tissue}/{pop}/{comp}.csv',tissue,pop,comp)
            opls = read_sig_dataset2(f'out/work/lipids/opls/{outlier}/{tissue}/{pop}/{comp}.csv',tissue,pop,comp)
            zscore = read_sig_dataset2(f'out/work/lipids/zscore/{outlier}/{tissue}/{pop}/{comp}.csv',tissue,pop,comp)
            glm['LMID'] = glm.apply(lambda u: get_lipidmap_id(u.name),axis=1)
            glm['Category'] = glm.apply(lambda u: ali.lipidmaps[u['LMID']]['CATEGORY'] if isinstance(u['LMID'],str) else nan,axis=1)
            glm['MainClass'] = glm.apply(lambda u: ali.lipidmaps[u['LMID']]['MAIN_CLASS'] if isinstance(u['LMID'],str) else nan,axis=1)
            glm['Saturation'] = glm.apply(lambda u: assignfa(u['LMID']) if isinstance(u['LMID'],str) else nan,axis=1)
            cross_cond_significance.append(glm)
cross_cond_significance = concat(cross_cond_significance,axis=0).dropna(subset=['LMID'])
cross_cond_significance = cross_cond_significance.rename({'Pr(>|z|)':'p-val','Estimate':'Slope'},axis=1)
cross_cond_significance.index.name = 'Name'
if args.out_cross_pop:
    cross_cond_significance.to_csv(args.out_cross_cond)
