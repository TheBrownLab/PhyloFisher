from peewee import *

database = SqliteDatabase(None)

class UnknownField(object):
    def __init__(self, *_, **__): pass

class BaseModel(Model):
    class Meta:
        database = database

class Genes(BaseModel):
    name = CharField(null=True, unique=True)

    class Meta:
        table_name = 'genes'

class Taxonomies(BaseModel):
    color = TextField(null=True)
    taxonomy = CharField(null=True, unique=True)

    class Meta:
        table_name = 'taxonomies'

class Metadata(BaseModel):
    data_type = TextField()
    higher_taxonomy = ForeignKeyField(column_name='higher_taxonomy_id', field='id', model=Taxonomies)
    long_name = TextField()
    lower_taxonomy = ForeignKeyField(backref='taxonomies_lower_taxonomy_set', column_name='lower_taxonomy_id', field='id', model=Taxonomies)
    short_name = CharField(unique=True)
    source = TextField()

    class Meta:
        table_name = 'metadata'

class Sequences(BaseModel):
    gene = ForeignKeyField(column_name='gene_id', field='id', model=Genes)
    header = TextField(index=True)
    is_paralog = BooleanField()
    metadata = ForeignKeyField(column_name='metadata_id', field='id', model=Metadata)
    sequence = TextField()

    class Meta:
        table_name = 'sequences'

class SqliteSequence(BaseModel):
    name = BareField(null=True)
    seq = BareField(null=True)

    class Meta:
        table_name = 'sqlite_sequence'
        primary_key = False

