from obj_tables import (Model, TableFormat, StringAttribute, UrlAttribute, OneToOneAttribute, ManyToOneAttribute, OneToManyAttribute, FloatAttribute, IntegerAttribute, BooleanAttribute)
import obj_tables
from obj_tables.io import Writer

class Address(obj_tables.Model):
    street = obj_tables.StringAttribute(unique=True, primary=True, verbose_name='Street')
    city = obj_tables.StringAttribute(verbose_name='City')
    state = obj_tables.StringAttribute(verbose_name='State')
    zip_code = obj_tables.StringAttribute(verbose_name='Zip code')
    country = obj_tables.StringAttribute(verbose_name='Country')

    class Meta(obj_tables.Model.Meta):
        table_format = obj_tables.TableFormat.multiple_cells
        attribute_order = ('street', 'city', 'state', 'zip_code', 'country',)
        verbose_name = 'Address'
        verbose_name_plural = 'Addresses'

class Company(obj_tables.Model):
    name = obj_tables.StringAttribute(unique=True, primary=True, verbose_name='Name')
    url = obj_tables.UrlAttribute(verbose_name='URL')
    address = obj_tables.OneToOneAttribute(Address, related_name='company', verbose_name='Address')

    class Meta(obj_tables.Model.Meta):
        table_format = obj_tables.TableFormat.column
        attribute_order = ('name', 'url', 'address',)
        verbose_name = 'Company'
        verbose_name_plural = 'Companies'

google = Company(name='Google',
                 url='https://www.google.com/',
                 address=Address(street='1600 Amphitheatre Pkwy',
                                 city='Mountain View',
                                 state='CA',
                                 zip_code='94043',
                                 country='US'))

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
