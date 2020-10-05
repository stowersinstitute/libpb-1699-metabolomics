from obj_tables import (Model, TableFormat, StringAttribute, UrlAttribute, OneToOneAttribute, ManyToOneAttribute, OneToManyAttribute, FloatAttribute, IntegerAttribute, BooleanAttribute)
from obj_tables.io import Writer

class Population(Model):
    name = StringAttribute(unique=True,primary=True,verbose_name='Name')

    class Meta(Model.Meta):
        table_format = TableFormat.column
        attribute_order = ('name',)
        verbose_name = 'Population'
        verbose_name_plural = 'Populations'

Pachon = Population(name='Pachon')
Tinaja = Population(name='Tinaja')
Surface = Population(name='Surface')

class Observation(Model):
    name = StringAttribute(unique=False,primary=False)
    population = ManyToOneAttribute(Population,related_name='observationsx',verbose_name='Population')

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('name', 'population',)
        verbose_name = 'Observation'
        verbose_name_plural = 'Observations'


obs = Observation(
  name = 'this',
  population = Pachon)
