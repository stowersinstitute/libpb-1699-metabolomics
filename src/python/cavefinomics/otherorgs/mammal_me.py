import attr
from pandas import read_csv, DataFrame

mammals = ['SugarGlider', 'Hedgehog', 'Horse', 'Pig', 'Cattle', 'Goat', 'Cat', 'Dog', 'Bear', 'Badger', 'Vervet', 'Macaque', 'Rabbit', 'Chipmunk', 'GuineaPig', 'DamaralandMoleRat', 'NakedMoleRat', 'WhiteFootedMouse', 'Hamster', 'Gerbil', 'SpinyMouse', 'Rat', 'HouseMouse']

organs = ['Brain', 'Heart', 'Kidney', 'Liver']

@attr.s(auto_attribs=True)
class MaMammalianMe:
    annotation_csv: str
    normalized_value_csv: str

    def __attrs_post_init__(self):
        self.value_table = read_csv(self.normalized_value_csv, index_col='Metabolite')
        self.annotation_table = read_csv(self.annotation_csv, usecols=['Metabolite', 'KEGG_id'], index_col='Metabolite')
        # entries with multiple KEGG ids can just use one
        def fix_kegg_id(c):
            try:
                if not ',' in c:
                    return c
                else:
                    return sorted(c.split(','))[0]
            except TypeError:
                pass # nans, ignore
        self.annotation_table['KEGG_id'] = self.annotation_table['KEGG_id'].apply(fix_kegg_id)

        self.value_table = self.value_table.join(self.annotation_table).dropna()
        #print(self.value_table)
