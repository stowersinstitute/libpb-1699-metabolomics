import attr
from pandas import read_csv, DataFrame, concat
import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
from pyrsistent import freeze, thaw, PMap, m, pmap, PVector, pvector
from pampy import match, _, HEAD, TAIL
from pprint import pprint
from json import load
from functools import reduce
from typing import Type, Text
from collections import defaultdict
from itertools import chain
import os.path

# from pampy import match, _


class MissingDataError(Exception):
    pass

#https://hmdb.ca/metabolites/HMDB0000011

extra_categories = {
    'Carbohydrates': [
        # Monosaccharides
        'C07326', # 1,5-Anhydroglucitol
        'C08352', # 6-deoxyglucose
        'C02273', # digalacturonic acid
        'C02336', # fructose
        'C01094', # fructose-1-phosphate
        'C05345', # fructose-6-phosphate
        # Polysaccharides
        'C03661', # 1-kestose
        'C00026', # alpha-ketoglutarate
        'C01904', # arabitol
        'C00492', # raffinose
        # carbohydrate-related
        'C00158', # citric acid
        'C00503', # erythritol
        'C00122', # fumaric acid
        'C01235', # galactinol
        'C01507', # galactitol
        'C01113', # galactose-6-phosphate
        'C00800', # gluconic acid
        'C00198', # gluconic acid lactone
        'C00221', # glucose
        'C00103', # glucose-1-phosphate
        'C00092', # glucose-6-phosphate
        'C00258', # glyceric acid
        'C01432', # lactic acid
        'C00532', # lyxitol
        'C00711', # malic acid
        'C01835', # maltotriose
        'C08243', # melezitose
        'C00074', # phosphoenolpyruvate
        'C00345', # phosphogluconic acid
        'C00474', # ribitol
        'C01685', # ribonic acid
        'C00117', # ribose-5-phosphate
        'C00818', # saccharic acid
        'C00597', # 3-phosphoglycerate
        ],
    'Lipids/fatty acids': [
        'C00989', # 4-hydroxybutyric acid
        'C02979', # beta-glycerolphosphate
        'C02277', # dodecanol
        'C05401', # glycerol-3-galactoside
        'C03189', # glycerol-alpha-phosphate
        'C00294', # inosine
        'C00130', # inosine-5'-monophosphate
        'C03546', # inositol-4-monophosphate
        'C19670', # oleamide
        'C16537', # pentadecanoic acid
        'C00022', # pyruvic acid
        ],
    'Nucleic acids': [
        'C00043', # UDP-N-acetylglucosamine
        'C00262', # hypoxanthine
        'C00295', # orotic acid
        'C00385', # xanthine
        ],
    'Peptides': [
        'C02721', # 2-aminobutyric acid
        'C05984', # 2-hydroxybutanoic acid
        'C02630', # 2-hydroxyglutaric acid
        'C00322', # 2-ketoadipic acid
        'C00233', # 2-ketoisocaproic acid
        'C01744', # 3-(4-hydroxyphenyl)propionic acid
        'C01042', # N-acetylaspartic acid
        'C02714', # N-acetylputrescine
        'C00438', # N-carbamoylaspartate
        'C05829', # N-carbamylglutamate
        'C00043', # UDP-N-acetylglucosamine
        'C00791', # creatinine
        'C00489', # glutaric acid
        'C02989', # methionine sulfoxide
        'C00547', # noradrenaline
        'C00138', # putrescine
        ],
    'Vitamins and cofactors': [
        'C00872', # aminomalonate
        'C12276', # pantothenic acid
        'C00376', # tocopherol alpha-
        ],
    'Misc.': [
        'C06104', # adipic acid
        'C00956', # alpha-aminoadipic acid
        'C06180', # anabasine
        'C08261', # azelaic acid
        'C03557', # ciliatine
        'C00160', # glycolic acid
        'C07272', # maleimide
        'C00802', # parabanic acid
        'C00009', # phosphate
        'C00013', # pyrophosphate
        'C00346', # phosphoethanolamine
        'C00493', # shikimic acid
        'C00232', # succinate semialdehyde
        'C00042', # succinic acid
        'C00059', # sulfuric acid
        'C00245', # taurine
        'C06337', # terephthalic acid
        'C16884', # threitol
        'C01620', # threonic acid
        'C02470', # xanthurenic acid
        # Thiol metabolism
        'C05422', # dehydroascorbic acid
        # Uric acid metabolites
        'C00156', # 4-hydroxybenzoate
        'C00167', # UDP-glucuronic acid
        'C00086', # urea
        'C00366', # uric acid
        # Aromatic metabolites
        'C02814', # Hydroxyhydroquinone
        'C02480', # 2,4-hexadienedioic acid, degradation product of benzeme (https://3dprint.nih.gov/discover/3dpx-004575)
        'C00180', # benzoic acid
        'C00530', # hydroquinone
        'C00146', # phenol
        'C06202', # salicylaldehyde
        ],
}


def merge(u,v):
    #return match ((isinstance(u,PMap),isinstance(v,PMap)),
                  #(True, True), lambda x: u.update_with(merge, v),
                  #(True, _), u,
                  #(_, True), v,
                  #_, None)
    r = match((u,v),
              (Text, Text), lambda x,y: u if u==v else pvector((u,v)),
              (Text, type(None)), u,
              (type(None), Text), v,
              _, None
              )
    if r is not None:
        return r
    r = match ((type(u),type(v)),
                  (Type[PVector], Type[PVector]), lambda x,y: u+v,
                  (Type[PVector], type(None)), u,
                  (type(None), Type[PVector]), v,
                  (Type[PMap], Type[PMap]), lambda x,y: u.update_with(merge, v),
                  (Type[PMap], _), u,
                  (_, Type[PMap]), v,
                  _, None)
    return r

def cutoff_at_depth(map,depth):
    if (depth > 1):
        return pmap({u:cutoff_at_depth(v,depth-1) for (u,v) in map.items()})
    else:
        return pmap({u:None for (u,v) in map.items()})

def include_only(map,keys):
    return pmap({u:v for (u,v) in map.items() if u in keys})

def get_in(map,keychain):
    if len(keychain) > 1:
        return get_in(map[keychain[0]], keychain[1:])
    else:
        return map[keychain[0]]

def has_in(map,keychain):
    if len(keychain) > 1:
        return keychain[0] in map and has_in(map[keychain[0]], keychain[1:])
    else:
        return keychain[0] in map

@attr.s(auto_attribs=True)
class AstyanaxMe:
    data_csv: str
    kegg_compounds_file: str
    hmdb_file: str
    sample_sheet_path: str

    def __attrs_post_init__(self):
        self.column_table = (
            read_csv(self.data_csv, skiprows=0, nrows=7)
            .iloc[:, 10:]
            .transpose()
            .reset_index()
        )
        self.row_table = read_csv(self.data_csv, skiprows=8, usecols=list(range(10)))
        self.row_table = self.row_table.rename(columns={'BinBase name':'Name'})
        self.data = read_csv(self.data_csv, skiprows=8).iloc[:, 11:]
        with open(self.kegg_compounds_file) as f:
            self.kegg_compounds = freeze(load(f))
        with open(self.hmdb_file) as f:
            hmdb = load(f)
        self.hmdb_category = {}
        for metabolite in hmdb:
            inchikey_stripped = metabolite['inchikey'].split('-')[0]
            if 'kegg_id' in metabolite:
                if 'sub_class' in metabolite:
                    self.hmdb_category[metabolite['kegg_id']] = metabolite['sub_class']
                elif 'super_class' in metabolite:
                    self.hmdb_category[metabolite['kegg_id']] = metabolite['super_class']
        #pprint(self.hmdb_category)


        colnames = self.column_table.iloc[0]
        self.column_table = self.column_table.iloc[1:, :]
        self.column_table.columns = colnames
        self.sample_sheet = read_csv(self.sample_sheet_path,index_col='Sample Label')
        sample_dict = self.sample_sheet['Mass (mg)'].to_dict()
        self.column_table['Mass (mg)'] = self.column_table['comment'].apply(lambda c: sample_dict[c] if c in sample_dict else np.nan)
        self.treatment_descriptors = []
        for d in self.column_table.iterrows():
            species = d[1]['species']
            if species == 'Pachon cavefish':
                species = 'Pachon'
            if species == 'Tinaja cavefish':
                species = 'Tinaja'
            if species == 'Astyanax mexicanus surface fish':
                species = 'Surface'

            organ = d[1]['organ']
            if organ == 'Muscles':
                organ = 'Muscle'

            treatment = d[1]['treatment'].split('-')[0].replace('Pachon','').replace('Tinaja','').replace('Surface','').strip()
            self.treatment_descriptors.append((species, organ, treatment))

        def get_kegg_pathway_cats():
            merged_pathways = reduce(merge, (self.kegg_compounds[u] for u in list(self.get_data_by_kegg_id().set_index('KEGG').index) if u in self.kegg_compounds))['PATHWAY']
            #pprint(thaw(merged_pathways))

            # chosen to balance categories
            pathway_priorities = {
                  'map00040': ("\nPentose/\nglucuronate",10), # Pentose / glucuronate interconversion
                  "map00230": ("Purine\nmetabolism",12),
                  "map00240": ("Pyrimidine\n",10), # Pyrimidine metabolism
                  'map01200': ('Carbon\nmetabolism',31),
                  'map01210': ('2-Oxocarboxylic\nacid',23),
                  'map01212': ('Fatty acids',9),
                  'map01220': ('Degradation\nof aromatics',18),
                  'map01230': ('Biosynth.\namino acids',37),
              }
            pathways_of_interest = include_only(
              merged_pathways, [u for u in pathway_priorities.keys()])
            #print(thaw(pathways_of_interest))

            def compounds_with_pathway(pathway):
                return [u for u in list(self.get_data_by_kegg_id().set_index('KEGG').index) if u in self.kegg_compounds and self.kegg_compounds[u] is not None and (not 'PATHWAY' in self.kegg_compounds[u] or pathway in self.kegg_compounds[u]['PATHWAY'])]

            compounds_by_pathway = {p:compounds_with_pathway(p) for p in pathways_of_interest}
            #pprint(compounds_by_pathway)

            pathways_by_compound = {
                c: list(sorted([u for u,v in compounds_by_pathway.items() if c in v],
                               key=lambda x: pathway_priorities[x][1],
                               reverse=False)) for c in list(self.get_data_by_kegg_id().set_index('KEGG').index)}
            # fill in at least one pathway for compounds that still have none
            if False:
                for cpd,pathways in pathways_by_compound.items():
                    if pathways is None or len(pathways) == 0:
                        if cpd in self.kegg_compounds and self.kegg_compounds[cpd] is not None and 'PATHWAY' in self.kegg_compounds[cpd]:
                            def a(h,t):
                                if pathways_by_compound[cpd] is None:
                                    pathways_by_compound[cpd] = [h]
                                else:
                                    pathways_by_compound[cpd] += [h]
                                return None
                            match(list(sorted(self.kegg_compounds[cpd]['PATHWAY'])),
                                  [HEAD, TAIL], a,
                                  _, None
                                  )
            # fill in missing values as "Misc"
            if True:
                for cpd,pathways in pathways_by_compound.items():
                    if pathways is None or len(pathways) == 0:
                        pathways_by_compound[cpd] = ['mapXXXXX']
            #pprint(pathways_by_compound)
            # missing carbon pathways
            pathways_by_compound['C00103'] += ['map00010', 'map01200'] # glucose-1-phosphate
            pathways_by_compound['C00221'] += ['map00010', 'map01200'] # beta D-glucose
            # missing fatty acid pathways
            pathways_by_compound['C00712'] += ['map01212'] # oleic acid
            #print(pathways_by_compound['C00103'])
            represented_pathways = list(sorted(set(chain.from_iterable(pathways_by_compound.values()))))
            #print(represented_pathways)
            pathways_by_compound_d = DataFrame([[1 if pathway in pathways else 0 for pathway in represented_pathways] for cpd,pathways in pathways_by_compound.items()], columns=represented_pathways, index=pathways_by_compound.keys())
            #print(pathways_by_compound_d.sum())
            return pathways_by_compound, pathway_priorities

        organize_by = 'hmdb'
        if organize_by == 'brite':
            merged_brite = reduce(merge, (self.kegg_compounds[u] for u in list(self.get_data_by_kegg_id().set_index('KEGG').index) if u in self.kegg_compounds))['BRITE']
            #pprint(thaw(cutoff_at_depth(merged_brite,1)))
            categories_of_interest = include_only(
                merged_brite, [
                    'Compounds with biological roles',
                    'Glycosides',
                    'Lipids',
                ])

            def compounds_with_brite_category(category):
                return [
                    u for u,v in self.kegg_compounds.items()
                    if v is not None and has_in(v, ['BRITE']+category)]

            self.compounds_by_category = {
                'Carbohydrates': compounds_with_brite_category(
                    ['Compounds with biological roles', 'Carbohydrates']),
                'Lipids/fatty acids': compounds_with_brite_category(
                    ['Compounds with biological roles', 'Lipids']),
                'Nucleic acids': compounds_with_brite_category(
                    ['Compounds with biological roles', 'Nucleic acids']),
                'Peptides': compounds_with_brite_category(
                    ['Compounds with biological roles', 'Peptides']),
                'Steroids': compounds_with_brite_category(
                    ['Compounds with biological roles', 'Steroids']),
                'Vitamins and cofactors': compounds_with_brite_category(
                    ['Compounds with biological roles', 'Vitamins and cofactors']),
            }

            for u,compounds in extra_categories.items():
                if u in self.compounds_by_category:
                    self.compounds_by_category[u] += compounds
                else:
                    self.compounds_by_category[u] = compounds

            self.classes_by_compound = {c: [u for u,v in self.compounds_by_category.items() if c in v] for c in list(self.get_data_by_kegg_id().set_index('KEGG').index)}

            self.compounds_by_category_from_dataset = defaultdict(list)
            for compound,v in self.classes_by_compound.items():
                for category in v:
                    self.compounds_by_category_from_dataset[category].append(compound)
        elif organize_by == 'pathway':
            pathways_by_compound, pathway_priorities = get_kegg_pathway_cats()
            self.compounds_by_category_from_dataset = defaultdict(list)
            for compound,pathways in pathways_by_compound.items():
                for pathway_id in pathways:
                    if pathway_id in pathway_priorities:
                        self.compounds_by_category_from_dataset[pathway_priorities[pathway_id][0]].append(compound)
                    else:
                        self.compounds_by_category_from_dataset['Misc.'].append(compound)
            #print(self.compounds_by_category_from_dataset)
        elif organize_by == 'hmdb':
            categories = defaultdict(list)
            for c in list(self.get_data_by_kegg_id().set_index('KEGG').index):
                if c in self.hmdb_category:
                    categories[self.hmdb_category[c]].append(c)
            #pprint(dict(categories))
            self.compounds_by_category_from_dataset = {}

            # manually assign nucleotides
            self.compounds_by_category_from_dataset['Nucleotides'] = \
              categories['Nucleosides, nucleotides, and analogues'] + \
              categories['Purine ribonucleotides'] + \
              categories['Purines and purine derivatives'] + categories['Pyridine carboxaldehydes'] + \
              categories['Pyridinecarboxylic acids and derivatives'] + \
              categories['Pyrimidine 2\'-deoxyribonucleosides'] + \
              categories['Pyrimidine nucleotide sugars'] + \
              categories['Pyrimidine ribonucleotides'] + \
              categories['Pyrimidines and pyrimidine derivatives']
            #nucleotides = set(self.compounds_by_category_from_dataset['Nucleotides'])
            # manually assign fatty acids
            self.compounds_by_category_from_dataset['Fatty\nacids'] = categories['Fatty acids and conjugates'] + categories['Fatty alcohols'] + categories['Fatty amides'] + categories['Lineolic acids and derivatives'] + categories[''] + [
              'C00489', # Glutarate, fatty acid degradation
              'C05984', # 2-Hydroxybutanoic acid, fatty acids in BRITE, insulin resistance
              'C08362', # palmitoleic acid
              'C00376', # retinal
            ]
            # carbon metabolism
            self.compounds_by_category_from_dataset['Carbohydrates /\nCCM'] = [
              'C00022', # Pyruvate
              'C00026', # 2-Oxoglutarate
              'C00042', # Succinate
              'C00074', # Phosphoenolpyruvate
              'C00122', # Fumarate
              'C00137', # myo-Inositol
              'C00158', # Citrate
              'C00199', # D-Ribulose 5-phosphate
              'C00800', # L-Gulonate
              'C03546', # myo-Inositol 4-phosphate
              'C00989', # 4-hydroxybutyric acid
              'C02336', # fructose
              'C01113', # galactose-6-phosphate
              'C00252', # isomaltose
              'C01432', # lactic acid
              'C00711', # malic acid
              'C08250', # sophorose
              'C00159', # mannose
              ] + \
                categories['Carbohydrates and carbohydrate conjugates']
            # amino acids
            self.compounds_by_category_from_dataset['Amino\nacids'] = \
              categories['Amino acids, peptides, and analogues'] + [
                'C00078', # Tryptophan
                'C00233', # 4-Methyl-2-oxopentanoate (valine / isoleucine degradation)
                'C00322', # 2-Oxoadipate
                'C02714', # N-Acetylputrescine
                'C00493', # Shikimate
                'C00547', # L-Noradrenaline
                'C02989', # methionine sulfoxide
              ]
            # misc. / secondary metabolites
            self.compounds_by_category_from_dataset['Misc. / sec.\nmetabolites'] = [
              'C00009', # Phosphate
              'C00013', # Diphosphate
              'C00059', # Sulfate
              'C00086', # Urea
              'C00156', # 4-Hydroxybenzoate
              'C00187', # Cholesterol
              'C00189', # Ethanolamine
              'C00245', # Taurine
              'C00346', # Ethanolamine phosphate
              'C00639', # Prostaglandin F2alpha, arachidonic acid metabolism
              'C01744', # Phloretate
              'C02470', # Xanthurenic acid
              'C02630', # 2-Hydroxyglutarate
              'C02979', # Glycerol 2-phosphate
              'C03557', # 2-Aminoethylphosphonate
              'C05332', # Phenethylamine
              'C05401', # Galactosylglycerol
              'C06180', # Anabasine
              'C06202', # Salicylaldehyde
              'C06337', # Terephthalate
              'C00530', # Hydroquinone
              'C03844', # Pinitol
              'C02814', # 1,2,4-benzenetriol
              'C02721', # 2-aminobutyric acid
              'C08352', # 6-deoxyglucose
              'C00072', # ascorbic acid
              'C05422', # dehydroascorbic acid
              'C00180', # benzoic acid
              'C00503', # erythritol
              'C00051', # glutathione
              'C03189', # glycerol-alpha-phosphate
              'C00160', # glycolic acid
              'C07272', # maleimide
              'C12276', # pantothenic acid
              'C00802', # parabanic acid
              'C00146', # phenol
              'C00138', #
              'C00134', # putrescine
            ]
            curated_categories = set(chain.from_iterable(self.compounds_by_category_from_dataset.values()))
            #print(curated_categories)


            # remove nucleotides from categories
            categories = {c:[cpd for cpd in cpds if cpd not in curated_categories] for c,cpds in categories.items()}
            # remove zero-length lists
            categories = {c:cpds for c,cpds in categories.items() if len(cpds) > 0}
            #pprint(categories)

            # split out low-membership categories
            short_categories = {c:cpds for c,cpds in categories.items() if len(cpds) <= 3}
            categories = {c:cpds for c,cpds in categories.items() if len(cpds) > 3}
            #pprint(short_categories)
            #categories['Amino acids, peptides, and analogues']
            recat_compounds = list(chain.from_iterable(cpds for c,cpds in short_categories.items()))
            #pprint(recat_compounds)
            pathways_by_compound, pathway_priorities = get_kegg_pathway_cats()
            #pprint({c:pathways_by_compound[c] for c in recat_compounds})


        self.kegg_to_name_map = self.get_kegg_to_name_map()
        self.kegg_to_name_map['C01080'] = 'Coenzyme F420'
        self.kegg_to_name_map['C0108'] = 'Coenzyme F420' # error in data - wrong num digits
        self.kegg_to_name_map['C00380'] = 'cytosine' # misspelled as 'cytosin'
        #pprint({f'{kegg_to_name_map[u]} ({u})': v for u,v in self.classes_by_compound.items()})

    def get_data_subset(
        self, *, kegg_metabolites, exp_property, property_values,
    ):
        rows = self.row_table["KEGG"].isin(kegg_metabolites)
        if not rows.any():
            return DataFrame()
        cols = self.column_table["treatment"].isin(property_values)

        result = self.data.loc[rows, list(cols)]
        result.index = kegg_metabolites
        result.columns = self.column_table["treatment"].loc[cols]
        return result

    def t_test_ind(self, *, kegg_metabolites, exp_property, property_values):
        from numpy import log
        import math

        def single_metabolite_test(kegg_compound):
            return ttest_ind(g[0][1][kegg_compound], g[1][1][kegg_compound]).pvalue

        d = self.get_data_subset(
            kegg_metabolites=kegg_metabolites,
            exp_property=exp_property,
            property_values=property_values,
        )
        if len(d.index) == 0:
            return DataFrame(data={c: math.inf for c in kegg_metabolites}, index=[0],)

        d = d.transpose().apply(log).groupby("treatment")
        if d.ngroups == 2:
            g = [g for g in d]
            return DataFrame(
                data={c: single_metabolite_test(c) for c in kegg_metabolites},
                index=[0],
            )
        else:
            raise RuntimeError(f"Too many groups for t-test: {d.ngroups}")

    def t_test_ind_combined(self, *, kegg_metabolites, exp_property, property_values):
        from numpy import log
        import scipy.stats

        d = self.t_test_ind(
            kegg_metabolites=kegg_metabolites,
            exp_property=exp_property,
            property_values=property_values,
        )
        ntests = len(kegg_metabolites)
        if ntests == 1:
            return d.iloc[0, 0]
        else:
            return scipy.stats.chi2.sf(
                -(2.0 * d * ntests / 80.0).apply(log).to_numpy().sum(), 2 * ntests
            )

    def get_twoway_log2_fc(
        self, *, kegg_metabolites, exp_property, property_values1, property_values2,
    ):
        """
        Returns the average (over all specified genes) of a two-way comparison between pop1 and pop2 for a given quantity.
        """
        from numpy import log2

        d1 = self.get_data_subset(
            kegg_metabolites=kegg_metabolites,
            exp_property=exp_property,
            property_values=property_values1,
        )
        stats1 = DataFrame({"Mean": d1.mean(axis=1), "Std": d1.std(axis=1)})

        d2 = self.get_data_subset(
            kegg_metabolites=kegg_metabolites,
            exp_property=exp_property,
            property_values=property_values2,
        )
        stats2 = DataFrame({"Mean": d2.mean(axis=1), "Std": d2.std(axis=1)})

        log2_mean_fc = DataFrame(
            {"log2fc": log2(stats2["Mean"]) - log2(stats1["Mean"])}
        )
        return float(log2_mean_fc.mean())

    def get_twoway_log_inverse_rel_cv(
        self, *, kegg_metabolites, exp_property, property_values1, property_values2,
    ):
        """
        Returns the average (over all specified genes) of a two-way comparison between pop1 and pop2 for a given quantity.
        """
        from numpy import log

        d1 = self.get_data_subset(
            kegg_metabolites=kegg_metabolites,
            exp_property=exp_property,
            property_values=property_values1,
        )
        stats1 = DataFrame({"Mean": d1.mean(axis=1), "Std": d1.std(axis=1)})

        d2 = self.get_data_subset(
            kegg_metabolites=kegg_metabolites,
            exp_property=exp_property,
            property_values=property_values2,
        )
        stats2 = DataFrame({"Mean": d2.mean(axis=1), "Std": d2.std(axis=1)})

        inverse_cvs = DataFrame(
            {
                "invCV1": stats1["Mean"] / stats1["Std"],
                "invCV2": stats2["Mean"] / stats2["Std"],
            }
        )

        log_inv_rel_cv = DataFrame(
            {"logInvRelCV": log(inverse_cvs["invCV2"]) - log(inverse_cvs["invCV1"])}
        )
        return float(log_inv_rel_cv.sum())

    def get_data_by_kegg_id(self):
        return concat((self.row_table['KEGG'], self.data), axis=1).dropna()

    def get_data_by_inchi_key(self):
        return concat((self.row_table['InChI Key'], self.data), axis=1).dropna()


    def get_data_by_name(self):
        return concat((self.row_table['Name'], self.data), axis=1).loc[self.row_table['KEGG'].notna(),:].dropna()


    def get_data_by_name_and_kegg(self):
        return concat((self.row_table['Name'], self.row_table['KEGG'], self.data), axis=1).loc[self.row_table['KEGG'].notna(),:].dropna()


    def get_name_kegg_columns(self):
        return concat((self.row_table['Name'], self.row_table['KEGG']), axis=1).dropna()


    def get_kegg_to_name_map(self):
        d = self.get_name_kegg_columns()
        return {str(kegg):str(name) for _, (name,kegg) in d.iterrows()}

    def get_name_to_kegg_map(self):
        d = self.get_name_kegg_columns()
        return {str(name):str(kegg) for _, (name,kegg) in d.iterrows()}
