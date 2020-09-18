from argparse import ArgumentParser
from cavefinomics import AstyanaxMe
import matplotlib.pyplot as plt
from matplotlib.ticker import EngFormatter
import matplotlib.ticker as plticker
from itertools import repeat, chain
from collections import defaultdict
flatten = chain.from_iterable
import os
import json
from pandas import concat, DataFrame, read_csv


parser = ArgumentParser(description="Plot metabolites of interest.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--exclude-outlier", type=bool, help="Exclude the outliers?")
parser.add_argument("--out-mtic", type=bool, help="Output for mTIC")
parser.add_argument("--out-cross-pop", type=bool, help="Output for cross-pop comparison")
parser.add_argument("--out-starvation-resp", type=bool, help="Output for starvation response")
args = parser.parse_args()

ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    sample_sheet_path=args.sample_sheet,
    hmdb_file=args.hmdb,
    )

kegg_to_hmdb = defaultdict(list)
with open(args.hmdb) as f:
    hmdb = json.load(f)
for metabolite in hmdb:
    #print(metabolite)
    if 'kegg_id' in metabolite:
        kegg_to_hmdb[metabolite['kegg_id']].append(metabolite['accession'])

if args.exclude_outlier:
    outlier = 'kein-Ausreißern'
else:
    outlier = 'mit-Ausreißern'

pops = ['Pachon', 'Tinaja', 'Surface']
tissues = ['Brain', 'Muscle', 'Liver']
conditions = ['4d Starved', '30d Starved', 'Refed']
conditions_really_short = ['4', r'30', 'R']
comparisons = {'PvS':('Pachon','Surface'), 'TvS':('Tinaja','Surface'), 'PvT':('Pachon','Tinaja')}
condmap = {'30d':'30d Starved', '4d':'4d Starved', 'Ref': 'Refed'}
categories = {"Aminoacids":'Amino acids',"Carbohydrates_-CCM": 'Carbohydrates / CCM',"Fattyacids":'FattRy acids',"Misc._-_sec.metabolites":'Misc',"Nucleotides":'Nucleotides'}

#datasets = []
#sig = {}
#up = {}
cross_pop_significance = []
def read_sig_dataset(filepath,cat,tissue,cond,comp):
    d = read_csv(filepath,index_col=0)
    d['Category'] = cat
    d['Tissue'] = tissue
    d['Condition'] = cond
    d['Comparison'] = comp
    return d

for cat in categories:
    for tissue in tissues:
        for cond in condmap:
            for comp in comparisons:
                glm = read_sig_dataset(f'out/work/primary/glm/singlefactor/{outlier}/{cat}/{tissue}/{cond}/{comp}.csv',cat,tissue,cond,comp)
                opls = read_sig_dataset(f'out/work/primary/opls/{outlier}/{cat}/{tissue}/{cond}/{comp}.csv',cat,tissue,cond,comp)
                zscore = read_sig_dataset(f'out/work/primary/zscore/{outlier}/{cat}/{tissue}/{cond}/{comp}.csv',cat,tissue,cond,comp)
                cross_pop_significance.append(glm)
cross_pop_significance = concat(cross_pop_significance,axis=0).dropna()
cross_pop_significance = cross_pop_significance.rename({'Pr(>|z|)':'p-val','Estimate':'Slope'},axis=1)
cross_pop_significance.index.name = 'Name'
#print(cross_pop_significance)
if args.out_cross_pop:
    cross_pop_significance.to_csv(args.out_cross_pop)

group_comparisons = ['30vR','4vR','30v4']
conserved_strv_resp_significance = []
for cat in categories:
    for tissue in tissues:
        for cond in condmap:
            for comp in group_comparisons:
                glm = read_sig_dataset(f'out/work/primary/glm/singlefactor/{outlier}/{cat}/{tissue}/CvS/{comp}.csv',cat,tissue,cond,comp)
                conserved_strv_resp_significance.append(glm)
conserved_strv_resp_significance = concat(conserved_strv_resp_significance,axis=0).dropna()
conserved_strv_resp_significance = conserved_strv_resp_significance.rename({'Pr(>|z|)':'p-val','Estimate':'Slope'})
conserved_strv_resp_significance.index.name = 'Name'
conserved_strv_resp_significance = conserved_strv_resp_significance.rename({'Pr(>|z|)':'p-val','Estimate':'Slope'},axis=1)
conserved_strv_resp_significance.index.name = 'Name'
print(conserved_strv_resp_significance)
#stop
if args.out_starvation_resp:
    conserved_strv_resp_significance.to_csv(args.out_starvation_resp)

astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
astyanax_data.columns = (' '.join((c,str(n))) for c,n in zip(astyanax_data.columns,chain.from_iterable(repeat(range(1,6+1),9*3))))
#astyanax_data = astyanax_data.rename(ame.get_kegg_to_name_map(), axis=0)

outliers = ['Tinaja Liver Refed 6', 'Pachon Muscle Refed 5', 'Pachon Liver 30d Starved 3']

#def process_outlier(exclude,subset):
    #for outlier in outliers:
        #if exclude and outlier in subset.columns:
            #subset = subset.loc[:,~subset.columns.str.contains(outlier)]
    #return subset

mtic_data = []
kegg_to_name_map = ame.get_kegg_to_name_map()
for pop in pops:
    for tissue in tissues:
        for condition in conditions:
                #sig_subset = cross_pop_significance.loc[
                  #(cross_pop_significance['Tissue'] == tissue) & (cross_pop_significance['Condition'] == condition)]
                subset = astyanax_data.loc[:,
                  astyanax_data.columns.str.contains(pop) & astyanax_data.columns.str.contains(tissue) & astyanax_data.columns.str.contains(condition)]
                #print(subset)
                for compound,row in subset.iterrows():
                    #print(str(compound),row)
                    #print(compound)
                    for observation,val in row.iteritems():
                        #print(observation)
                        is_outlier = observation in outliers
                        rep = observation.split(' ')[-1]
                        mtic_data.append({
                          'Name': kegg_to_name_map[compound],
                          'KEGG': compound,
                          'HMDB': ','.join(kegg_to_hmdb[compound]) if compound in kegg_to_hmdb and len(kegg_to_hmdb[compound]) > 0 else None,
                          'Population': pop,
                          'Tissue': tissue,
                          'Condition': condition,
                          'Replicate': rep,
                          'Outlier': is_outlier,
                          'Raw_mTIC': val,
                        })

#mtic_data = concat(mtic_data,axis=0)
mtic_data = DataFrame(mtic_data).set_index('Name')
#print(mtic_data)
if args.out_mtic:
    conserved_strv_resp_significance.to_csv(args.out_mtic)
