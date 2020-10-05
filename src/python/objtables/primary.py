from pandas import read_csv
from obj_tables import (Model, TableFormat, StringAttribute, UrlAttribute, OneToOneAttribute, ManyToOneAttribute, OneToManyAttribute, FloatAttribute, IntegerAttribute, BooleanAttribute)
from obj_tables.io import Writer
from pprint import pprint

# because Pyhthon's defaultdict is worthless
#https://stackoverflow.com/questions/19399032/accessing-key-in-factory-of-defaultdict
class RealDefaultDict(dict):
    def __init__(self, factory):
        self.factory = factory
    def __missing__(self, key):
        self[key] = self.factory(key)
        return self[key]

#https://sandbox.karrlab.org/notebooks/obj_tables/1.%20Building%20and%20visualizing%20schemas.ipynb
#https://sandbox.karrlab.org/notebooks/obj_tables/2.%20Building%2C%20querying%2C%20editing%2C%20comparing%2C%20normalizing%2C%20and%20validating%20datasets.ipynb
#https://sandbox.karrlab.org/notebooks/obj_tables/3.%20Importing%2C%20exporting%2C%20converting%2C%20and%20pretty%20printing%20datasets.ipynb

data = read_csv('out/work/primary/merged-mtic.csv', index_col=0)
print(data)

class Population(Model):
    name = StringAttribute(unique=True,primary=True,verbose_name='Name')

    class Meta(Model.Meta):
        table_format = TableFormat.column
        attribute_order = ('name',)
        verbose_name = 'Population'
        verbose_name_plural = 'Populations'

class Tissue(Model):
    name = StringAttribute(unique=True,primary=True)

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('name',)
        verbose_name = 'Tissue'
        verbose_name_plural = 'Tissues'

class Condition(Model):
    name = StringAttribute(unique=True,primary=True)

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('name',)
        verbose_name = 'Condition'
        verbose_name_plural = 'Conditions'

class HMDBRecord(Model):
    hmdb_id = StringAttribute(unique=True,primary=True)

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('hmdb_id',)
        verbose_name = 'HMDBRecord'
        verbose_name_plural = 'HMDBRecords'

class ChEBIRecord(Model):
    chebi_id = StringAttribute(unique=True,primary=True)

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('chebi_id',)
        verbose_name = 'ChEBIRecord'
        verbose_name_plural = 'ChEBIRecords'

class Observation(Model):
    name = StringAttribute(unique=False,primary=False)
    kegg = StringAttribute(unique=False,primary=False)
    hmdb = OneToManyAttribute(HMDBRecord,related_name='observations')
    chebi = OneToManyAttribute(ChEBIRecord,related_name='observations')
    category = StringAttribute(unique=False,primary=False)
    population = ManyToOneAttribute(Population,related_name='observationsx',verbose_name='Population')
    tissue = ManyToOneAttribute(Tissue,related_name='observations')
    condition = ManyToOneAttribute(Condition,related_name='observations')
    replicate = IntegerAttribute()
    outlier = BooleanAttribute()
    raw_mtic = FloatAttribute()

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('name', 'kegg', 'hmdb', 'chebi', 'category', 'population', 'tissue', 'condition', 'replicate', 'outlier', 'raw_mtic',)
        verbose_name = 'Observation'
        verbose_name_plural = 'Observations'

def parse_list(u):
    if isinstance(u,str):
        return u.split(';')
    else:
        return []

Pachon = Population(name='Pachon')
Tinaja = Population(name='Tinaja')
Surface = Population(name='Surface')

Brain = Tissue(name='Brain')
Muscle = Tissue(name='Muscle')
Liver = Tissue(name='Liver')

cond_4d = Condition(name='4d Starved')
cond_30d = Condition(name='30d Starved')
cond_ref = Condition(name='Refed')

observations = []
hmdb = RealDefaultDict(lambda hmdb_id: HMDBRecord(hmdb_id=hmdb_id))
chebi = RealDefaultDict(lambda chebi_id: ChEBIRecord(chebi_id=chebi_id))
for i,row in data.iloc[:10,:].iterrows():
    hmdbs = {u: hmdb[u] for u in parse_list(row["HMDB"])}
    hmdb.update(hmdbs)

    chebis = {u: chebi[u] for u in parse_list(row["ChEBI"])}
    chebi.update(chebis)

    #if row['Population'] == 'Pachon':
        #pop = Pachon
    #elif row['Population'] == 'Tinaja':
        #pop = Tinaja
    #elif row['Population'] == 'Surface':
        #pop = Surface
    #else:
        #assert False

    if row['Tissue'] == 'Brain':
        tissue = Brain
    elif row['Tissue'] == 'Muscle':
        tissue = Muscle
    elif row['Tissue'] == 'Liver':
        tissue = Liver
    else:
        assert False

    if row['Condition'] == '4d Starved':
        condition = cond_4d
    elif row['Condition'] == '30d Starved':
        condition = cond_30d
    elif row['Condition'] == 'Refed':
        condition = cond_ref
    else:
        assert False

    obs = Observation(
      name = i.strip(),
      kegg = row['KEGG'],
      hmdb = hmdbs.values(),
      chebi = chebis.values(),
      category = row['Category'],
      population = Pachon,
      #tissue = tissue,
      #condition = condition,
      replicate = row['Replicate'],
      outlier = row['Outlier'],
      raw_mtic = row['Raw_mTIC'])

    observations.append(obs)

#schema = type('primary', (types.ModuleType, ), {
    #'Observation': Observation,
    #'models': [Observation],
#})

Writer().run('/tmp/x/*.csv', observations, models=[Population,Tissue,Condition,HMDBRecord,ChEBIRecord,Observation], write_toc=True, write_schema=True)
