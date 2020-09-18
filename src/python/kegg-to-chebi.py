from argparse import ArgumentParser
from cavefinomics import AstyanaxMe
from swagger_codegen.api.adapter.requests import RequestsAdapter
from cts import new_client, Configuration
from itertools import repeat, chain
flatten = chain.from_iterable
import json

parser = ArgumentParser(description="Plot metabolites of interest.")
parser.add_argument("--astyanax", type=str, help="Astyanax metabolomics csv file.")
parser.add_argument("--compounds", type=str, help="KEGG compounds file.")
parser.add_argument("--sample-sheet", type=str, help="Sample sheet.")
parser.add_argument("--hmdb", type=str, help="HMDB file.")
parser.add_argument("--out", type=str, help="Output")
args = parser.parse_args()

ame = AstyanaxMe(
    data_csv=args.astyanax,
    kegg_compounds_file=args.compounds,
    sample_sheet_path=args.sample_sheet,
    hmdb_file=args.hmdb,
    )

astyanax_data = ame.get_data_by_kegg_id().set_index('KEGG')
astyanax_data.columns = [' '.join(u) for u in ame.treatment_descriptors]
astyanax_data = astyanax_data.loc[:,['pools' not in c for c in astyanax_data.columns]]
astyanax_data.columns = (' '.join((c,str(n))) for c,n in zip(astyanax_data.columns,flatten(repeat(range(1,6+1),9*3))))

client = new_client(RequestsAdapter(), Configuration(host="https://cts.fiehnlab.ucdavis.edu"))

kegg_to_chebi = {}
keggs = list(sorted(astyanax_data.index))
for kegg in keggs:
    try:
        chebis = client.ids.convert("KEGG","ChEBI",kegg)[0].results
        kegg_to_chebi[kegg] = chebis
    except:
        pass

with open(args.out,'w') as f:
    json.dump(kegg_to_chebi,f)
