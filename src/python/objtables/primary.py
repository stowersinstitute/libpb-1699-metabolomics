from pandas import read_csv
from obj_tables import (Model, TableFormat, StringAttribute, UrlAttribute, OneToOneAttribute, ManyToOneAttribute, OneToManyAttribute, FloatAttribute)
from obj_tables.io import Writer
from pprint import pprint

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

class HMDBRecord(Model):
    hmdb_id = StringAttribute(unique=True,primary=True)

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('hmdb_id',)
        verbose_name = 'HMDBRecord'
        verbose_name_plural = 'HMDBRecords'

class Observation(Model):
    name = StringAttribute(unique=False,primary=False)
    kegg = StringAttribute(unique=False,primary=False)
    hmdb = OneToManyAttribute(HMDBRecord,related_name='observations')
    raw_mtic = FloatAttribute()

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('name', 'kegg', 'hmdb', 'raw_mtic',)
        verbose_name = 'Observation'
        verbose_name_plural = 'Observations'

def parse_list(u):
    if isinstance(u,str):
        return u.split(';')
    else:
        return []

observations = []
#hmdb = RealDefaultDict(lambda hmdb_id: HMDBRecord(hmdb_id))
hmdb = {}
for i,row in data.iloc[:10,:].iterrows():
    #print(row)
    #print(list(parse_list(row["HMDB"])))
    #hmdbs = [hmdb[u] for u in parse_list(row["HMDB"])]
    hmdbs = {u: hmdb[u] if u in hmdb else HMDBRecord(hmdb_id=u) for u in parse_list(row["HMDB"])}
    hmdb.update(hmdbs)
    print(hmdbs)
    obs = Observation(
      name = i.strip(),
      kegg = row['KEGG'],
      hmdb = hmdbs.values(),
      raw_mtic = row['Raw_mTIC'])
    #pprint(obs)
    observations.append(obs)

#schema = type('primary', (types.ModuleType, ), {
    #'Observation': Observation,
    #'models': [Observation],
#})

Writer().run('/tmp/x/*.csv', observations, models=[HMDBRecord,Observation], write_toc=True, write_schema=True)
