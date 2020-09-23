from argparse import ArgumentParser
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
parser.add_argument("--out", type=str, help="Output")
args = parser.parse_args()

with open(args.lipidmaps_fa) as f:
  lmfa = json.load(f)

with open(args.lipidmaps_json) as f:
    lm = json.load(f)

ali = AstyanaxLi(
    lipids_normalized=args.lipids_normalized,
    lipidmaps_js=args.lipidmaps_json,
    )

astyanax_data = ali.lmdata.iloc[:,ali.non_numeric_cols-1:]
astyanax_data.columns = (c[0] + ' ' + c[-2] + ' ' + ' '.join(c[1:-2]) + ' ' + c[-1] for c in (cc.split() for cc in astyanax_data.columns))
astyanax_data = astyanax_data.reindex(sorted(astyanax_data.columns), axis=1)

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
polarities = ['positive','negative']
conditions = ['4d Starved', '30d Starved', 'Refed']
conditions_really_short = ['4', r'30', 'R']
comparisons = {'PvS':('Pachon','Surface'), 'TvS':('Tinaja','Surface'), 'PvT':('Pachon','Tinaja')}
condmap = {'30d':'30d Starved', '4d':'4d Starved', 'Ref': 'Refed'}
categories = {"Aminoacids":'Amino acids',"Carbohydrates_-CCM": 'Carbohydrates / CCM',"Fattyacids":'FattRy acids',"Misc._-_sec.metabolites":'Misc',"Nucleotides":'Nucleotides'}

outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def get_lipidmap_id(inchikey):
    if inchikey in ali.lipidmaps_inchikey:
        return ali.lipidmaps_inchikey[inchikey]
    elif '-'.join(inchikey.split('-')[:2]) in ali.lipidmaps_inchikey2:
        return ali.lipidmaps_inchikey2['-'.join(inchikey.split('-')[:2])]
    elif inchikey.split('-')[0] in ali.lipidmaps_inchikey1:
        return ali.lipidmaps_inchikey1[inchikey.split('-')[0]]
    return nan

def fix_cols(df,tissue):
    d = df.copy()
    def j(c):
        return ' '.join(c)
    d.columns = (f'{c[0]} {tissue} {j(c[1:])} {n}' for c,n in zip((cc.split() for cc in d.columns),list(range(1,7))*27))
    d = process_outlier(d)
    d = d.reindex(sorted(d.columns), axis=1)
    return d

outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

def process_outlier(subset):
    for o in outliers:
        subset.loc[:,subset.columns.str.contains(o)] = nan
    return subset

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

def get_lipidmap_name(lmid):
    if 'NAME' in lm[lmid]:
        return lm[lmid]['NAME']
    else:
        return nan

# raw mtic data
mtic_data= []
for tissue in tissues:
    for polarity in polarities:
        d = ali.normalized[tissue,polarity]
        d = fix_cols(d,tissue)
        for pop in pops:
            for condition in conditions:
                subset = d.loc[:,
                  d.columns.str.contains(pop) & d.columns.str.contains(condition)]
                for inchikey,row in subset.iterrows():
                    for observation,val in row.iteritems():
                        is_outlier = observation in outliers
                        rep = observation.split(' ')[-1]
                        lmid = get_lipidmap_id(inchikey)
                        if isinstance(lmid,str):
                            name = get_lipidmap_name(lmid)
                            mtic_data.append({
                              'LMID': lmid,
                              'Name': name,
                              'InChIKey': inchikey,
                              'Category': ali.lipidmaps[lmid]['CATEGORY'] if isinstance(lmid,str) else nan,
                              'MainClass': ali.lipidmaps[lmid]['MAIN_CLASS'] if isinstance(lmid,str) else nan,
                              'Saturation': assignfa(lmid) if isinstance(lmid,str) else nan,
                              'Population': pop,
                              'Tissue': tissue,
                              'Condition': condition,
                              'Polarity': polarity[:1].capitalize()+polarity[1:],
                              'Replicate': rep,
                              'Outlier': is_outlier,
                              'Raw_mTIC': val,
                            })

mtic_data = DataFrame(mtic_data).set_index('LMID')
if args.out:
    mtic_data.to_csv(args.out)
