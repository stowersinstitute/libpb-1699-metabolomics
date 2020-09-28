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

parser = ArgumentParser(description="Lipids merged.")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument("--lipidmaps-json", type=str, help="Lipidmaps JSON.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--out-cross-pop", type=str, help="Output (cross-pop)")
args = parser.parse_args()

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

cross_pop_significance = []
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

# cross population comparison
for tissue in tissues:
    for cond in condmap:
        for comp in comparisons:
            glm = read_sig_dataset(f'out/work/lipids/glm/singlefactor/{outlier}/{tissue}/{cond}/{comp}.csv',tissue,cond,comp)
            opls = read_sig_dataset(f'out/work/lipids/opls/{outlier}/{tissue}/{cond}/{comp}.csv',tissue,cond,comp)
            zscore = read_sig_dataset(f'out/work/lipids/zscore/{outlier}/{tissue}/{cond}/{comp}.csv',tissue,cond,comp)
            glm['LMID'] = glm.apply(lambda u: get_lipidmap_id(u.name),axis=1)
            cross_pop_significance.append(glm)
cross_pop_significance = concat(cross_pop_significance,axis=0).dropna()
cross_pop_significance = cross_pop_significance.rename({'Pr(>|z|)':'p-val','Estimate':'Slope'},axis=1)
cross_pop_significance.index.name = 'Name'
if args.out_cross_pop:
    cross_pop_significance.to_csv(args.out_cross_pop)
