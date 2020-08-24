from argparse import ArgumentParser
import os

from numpy import array, log10
from pandas import DataFrame, Series
from sklearn.preprocessing import StandardScaler, Normalizer
import matplotlib.pyplot as plt
import matplotlib
from seaborn import heatmap
from sklearn.decomposition import PCA as PCA
from functools import reduce
from scipy.stats import chisquare, fisher_exact
import math
from pprint import pprint
from collections import defaultdict

from cavefinomics import AstyanaxMe, AstyanaxLi

parser = ArgumentParser(description="PCA for primary.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--lipids-normalized", type=str, help="Normalized lipids dir.")
parser.add_argument("--lipidmaps-json", type=str, help="Lipidmaps JSON.")
args = parser.parse_args()

ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    hmdb_file=args.hmdb,
    )

ali = AstyanaxLi(
    lipids_normalized=args.lipids_normalized,
    lipidmaps_js=args.lipidmaps_json,
    )


# get data and compute common metabolites
primary_data = ame.get_data_by_kegg_id().set_index('KEGG')
lipid_data = ali.lmdata.iloc[:,ali.non_numeric_cols-1:]
#print(astyanax_data)
print(f'{len(primary_data.index)} metabolites linked to KEGG')
print(f'{len(lipid_data.index)} lipids linked to LIPIDMAPS')
