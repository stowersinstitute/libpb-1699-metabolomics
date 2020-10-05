from obj_tables import (Model, TableFormat,
                        StringAttribute, ManyToOneAttribute)

class Company(Model):
    name = StringAttribute(unique=True, primary=True, verbose_name='Name')

    class Meta(Model.Meta):
        table_format = TableFormat.column
        attribute_order = ('name',)
        verbose_name = 'Company'
        verbose_name_plural = 'Companies'

class Person(Model):
    name = StringAttribute(unique=True, primary=True, verbose_name='Name')
    company = ManyToOneAttribute(Company, related_name='employees', verbose_name='Company')

    class Meta(Model.Meta):
        table_format = TableFormat.row
        attribute_order = ('name', 'company')
        verbose_name = 'Person'
        verbose_name_plural = 'People'


google = Company(name='Google')

pichai = Person(name='Sundar Pichai',
                company=google)
