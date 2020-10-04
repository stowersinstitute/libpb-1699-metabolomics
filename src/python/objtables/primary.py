from pandas import read_csv
from obj_tables import (Model, TableFormat, StringAttribute, UrlAttribute, OneToOneAttribute, ManyToOneAttribute, FloatAttribute)
from obj_tables.io import Writer
from pprint import pprint

#https://sandbox.karrlab.org/notebooks/obj_tables/1.%20Building%20and%20visualizing%20schemas.ipynb
#https://sandbox.karrlab.org/notebooks/obj_tables/2.%20Building%2C%20querying%2C%20editing%2C%20comparing%2C%20normalizing%2C%20and%20validating%20datasets.ipynb
#https://sandbox.karrlab.org/notebooks/obj_tables/3.%20Importing%2C%20exporting%2C%20converting%2C%20and%20pretty%20printing%20datasets.ipynb

data = read_csv('out/work/primary/merged-mtic.csv', index_col=0)
print(data)

class Observation(Model):
    name = StringAttribute(unique=False,primary=False)
    kegg = StringAttribute(unique=False,primary=False)
    raw_mtic = FloatAttribute()

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('name', 'kegg', 'raw_mtic',)
        verbose_name = 'Observation'
        verbose_name_plural = 'Observations'

def parse_list(u):
    if isinstance(u,str):
        return u.split(';')
    else:
        return []

observations = []
for i,row in data.iloc[:10,:].iterrows():
    #print(row)
    obs = Observation(
      name = i,
      kegg = row['KEGG'],
      raw_mtic = row['Raw_mTIC'])
    #pprint(obs)
    observations.append(obs)

#schema = type('primary', (types.ModuleType, ), {
    #'Observation': Observation,
    #'models': [Observation],
#})

Writer().run('/tmp/*.csv', observations, models=[Observation], write_toc=True, write_schema=True)
