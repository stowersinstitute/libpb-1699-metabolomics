from obj_tables import (Model, TableFormat, StringAttribute, UrlAttribute, OneToOneAttribute, ManyToOneAttribute, OneToManyAttribute, FloatAttribute, IntegerAttribute, BooleanAttribute)
import obj_tables
from obj_tables.io import Writer

class Company(obj_tables.Model):
    name = obj_tables.StringAttribute(unique=True, primary=True, verbose_name='Name')

    class Meta(Model.Meta):
        table_format = obj_tables.TableFormat.column
        attribute_order = ('name',)
        verbose_name = 'Company'
        verbose_name_plural = 'Companies'

google = Company(name='Google')

class Observation(Model):
    name = StringAttribute(unique=True,primary=True)
    company = ManyToOneAttribute(Company,related_name='observationsx',verbose_name='Population')

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('name', 'company',)
        verbose_name = 'Observation'
        verbose_name_plural = 'Observations'


obs = Observation(
  name = 'this',
  company = google)
