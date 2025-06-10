#!/usr/bin/env python
import configparser
import os
import sys
import textwrap
from pathlib import Path
import pandas as pd
from peewee import *
from phylofisher import help_formatter
from phylofisher.utilities import build_database
from phylofisher.db_map import database, Taxonomies, Metadata, Sequences

pd.options.display.float_format = '{:,.0f}'.format


class MetadataHandler:
    '''
    A class to query metadata, orthologs, and paralogs from SQLite database via Peewee ORM.
    '''
    
    def __init__(self):
        '''
        Initialize the MetadataHandler by querying the Metadata table
        '''
        higher = Taxonomies.alias('higher')
        lower = Taxonomies.alias('lower')

        query = (
            Metadata
            .select(
                Metadata.short_name.alias('Unique ID'),
                Metadata.long_name.alias('Long Name'),
                Metadata.data_type.alias('Data Type'),
                Metadata.source.alias('Source'),
                higher.taxonomy.alias('Higher Taxonomy'),
                lower.taxonomy.alias('Lower Taxonomy'),
            )
            .join(higher, on=(Metadata.higher_taxonomy == higher.id))
            .switch(Metadata)
            .join(lower, on=(Metadata.lower_taxonomy == lower.id))
        )

        self.metadata = pd.DataFrame(list(query.dicts()))
        self.metadata['Orthologs'] = self.metadata['Unique ID'].map(self.count_sequences(is_paralog=False))
        self.metadata['Paralogs'] = self.metadata['Unique ID'].map(self.count_sequences(is_paralog=True))

        self.org_meta = self.metadata[['Unique ID', 'Long Name', 'Higher Taxonomy', 'Lower Taxonomy',
                                    'Data Type', 'Orthologs', 'Paralogs', 'Source']]
        self.metadata = self.metadata[['Unique ID', 'Long Name', 'Higher Taxonomy', 'Lower Taxonomy',
                                    'Data Type', 'Orthologs', 'Paralogs']]

    def count_sequences(self, is_paralog=False):
        '''
        Count the number of orthologs or paralogs for each organism in the database.

        :param is_paralog: If True, count paralogs; if False, count orthologs; default is False.
        :type is_paralog: bool, optional
        :return: A dictionary mapping organism short names to the count of orthologs or paralogs.
        :rtype: dict
        '''
        q = (
            Sequences
            .select(
                Metadata.short_name.alias('uid'),
                fn.COUNT(Sequences.id).alias('count')
            )
            .join(Metadata)
            .where(Sequences.is_paralog == is_paralog)
            .group_by(Metadata.short_name)
        )
        return {row['uid']: row['count'] for row in q.dicts()}

    def higher_taxonomy(self):
        '''
        Group metadata by higher taxonomy and count the number of organisms in each group.

        :return: A DataFrame with higher taxonomy and the count of organisms.
        :rtype: pandas.DataFrame
        '''
        return self.metadata.groupby(['Higher Taxonomy']
                                     ).size().reset_index().rename(columns={0: 'Organisms'})

    def lower_taxonomy(self):
        '''
        Group metadata by lower taxonomy and count the number of organisms in each group.

        :return: A DataFrame with lower taxonomy and the count of organisms.
        :rtype: pandas.DataFrame
        '''
        return self.metadata.groupby(['Higher Taxonomy', 'Lower Taxonomy']
                                     ).size().reset_index().rename(columns={0: 'Organisms'})

    def get_higher(self, term):
        '''
        Get metadata for a specific higher taxonomy.

        :param term: The higher taxonomy to filter by.
        :type term: str
        :return: A DataFrame containing metadata for the specified higher taxonomy.
        :rtype: pandas.DataFrame
        '''
        try:
            ret = self.metadata[self.metadata['Higher Taxonomy'] == term].sort_values(['Lower Taxonomy', 'Unique ID']).reset_index(drop=True)
        except KeyError:
            sys.exit(f'{term} is not a valid Higher Taxonomy')
        return ret

    def get_lower(self, term):
        '''
        Get metadata for a specific lower taxonomy.

        :param term: The lower taxonomy to filter by.
        :type term: str
        :return: A DataFrame containing metadata for the specified lower taxonomy.
        :rtype: pandas.DataFrame
        '''
        try:
            ret = self.metadata[self.metadata['Lower Taxonomy'] == term].sort_values('Unique ID').reset_index(drop=True)
        except KeyError:
            sys.exit(f'{term} is not a valid Lower Taxonomy')
        return ret

    def get_org(self, term):
        '''
        Get metadata for a specific organism by its unique ID.

        :param term: The unique ID of the organism to filter by.
        :type term: str
        :return: A DataFrame containing metadata for the specified organism.
        :rtype: pandas.DataFrame
        '''
        org_df = self.org_meta.copy().set_index('Unique ID')
        try:
            df = pd.DataFrame([org_df.loc[term, :]]).transpose()
        except KeyError:
            sys.exit(f'{term} is not a valid UniqueID')
        return df

def get_or_create_taxonomy(name):
    '''
    Get or create a Taxonomy entry in the database.

    :param name: The name of the taxonomy to get or create.
    :type name: str
    :return: The Taxonomy object corresponding to the given name.
    :rtype: Taxonomies
    '''
    return Taxonomies.get_or_create(taxonomy=name, defaults={"color": None})[0]

def update_metadata_table(tsv_path, dry_run=False):
    '''
    Update the Metadata table in the database using a TSV file.

    :param tsv_path: Path to the TSV file containing metadata updates.
    :type tsv_path: str
    :param dry_run: If True, do not apply changes to the database, just print what would be changed; defaults to False.
    :type dry_run: bool, optional
    '''

    # Read the TSV file into a DataFrame
    df = pd.read_csv(tsv_path, sep="\t", header=0)
    
    # Check if the required columns are present in the DataFrame
    required_columns = ['Unique ID', 'Long Name', 'Higher Taxonomy', 'Lower Taxonomy', 'Data Type', 'Source']
    if not all(col in df.columns for col in required_columns):
        sys.exit(f"Error: TSV file must contain the following columns: {', '.join(required_columns)}")

    # Loop through each row in the DataFrame and update the Metadata table
    for _, row in df.iterrows():
        try:
            meta = Metadata.get(Metadata.short_name == row['Unique ID'])
            ht = get_or_create_taxonomy(row['Higher Taxonomy'])
            lt = get_or_create_taxonomy(row['Lower Taxonomy'])

            changes = {}

            if meta.long_name != row['Long Name']:
                changes['long_name'] = (meta.long_name, row['Long Name'])

            if meta.data_type != row['Data Type']:
                changes['data_type'] = (meta.data_type, row['Data Type'])

            if meta.source != row['Source']:
                changes['source'] = (meta.source, row['Source'])

            if meta.higher_taxonomy.id != ht.id:
                changes['higher_taxonomy'] = (meta.higher_taxonomy.taxonomy, ht.taxonomy)

            if meta.lower_taxonomy.id != lt.id:
                changes['lower_taxonomy'] = (meta.lower_taxonomy.taxonomy, lt.taxonomy)

            if changes:
                if dry_run:
                    print(f"\n[Dry Run] Changes for {meta.short_name}:")
                    for field, (old, new) in changes.items():
                        print(f"  {field}: '{old}' -> '{new}'")
                else:
                    meta.long_name = row['Long Name']
                    meta.data_type = row['Data Type']
                    meta.source = row['Source']
                    meta.higher_taxonomy = ht
                    meta.lower_taxonomy = lt
                    meta.save()
                    print(f"Metadata table updated using {tsv_path}")
            else:
                print(f"Skipping {row['Unique ID']}: No changes detected.")


        except Metadata.DoesNotExist:
            print(f"Skipping: {row['Unique ID']} not found in metadata table.")
        except Taxonomies.DoesNotExist as e:
            print(f"Skipping {row['Unique ID']}: taxonomy not found - {e}")

def show_taxonomy_colors():
    '''
    Print all taxonomy terms and their associated colors.

    :return: None
    '''
    rows = Taxonomies.select().order_by(Taxonomies.taxonomy.asc())
    print(f"{'Taxonomy':<30} {'Color'}")
    for row in rows:
        color = row.color if row.color is not None else "NULL"
        print(f"{row.taxonomy:<30} {color}")


def update_taxonomy_colors(tsv_path, dry_run=False):
    '''
    Update the color field of taxonomies in the database from a TSV file.

    :param tsv_path: Path to the TSV file with columns Taxonomy and Color.
    :type tsv_path: str
    :param dry_run: If True, print changes but do not apply them to the database.
    :type dry_run: bool
    '''
    df = pd.read_csv(tsv_path, sep="\t")
    required_columns = ['Taxonomy', 'Color']
    if not all(col in df.columns for col in required_columns):
        sys.exit(f"Error: TSV file must contain the following columns: {', '.join(required_columns)}")

    for _, row in df.iterrows():
        try:
            tax = Taxonomies.get(Taxonomies.taxonomy == row['Taxonomy'])
            new_color = row['Color'] if pd.notnull(row['Color']) else None

            if tax.color != new_color:
                if dry_run:
                    print(f"[Dry Run] {tax.taxonomy}: '{tax.color}' -> '{new_color}'")
                else:
                    tax.color = new_color
                    tax.save()
                    print(f"Updated {tax.taxonomy}: '{tax.color}'")
            else:
                print(f"Skipping {tax.taxonomy}: Color is already '{tax.color}'.")

        except Taxonomies.DoesNotExist:
            print(f"Skipping: Taxonomy '{row['Taxonomy']}' not found in database.")


def update_unique_ids(threads, tsv_path, dry_run=False):
    '''
    Update organism Unique IDs in the Metadata table and reflect changes in Sequences headers.

    :param tsv_path: Path to the TSV file with columns Old ID and New ID.
    :type tsv_path: str
    :param dry_run: If True, just print the planned changes.
    :type dry_run: bool
    '''
    df = pd.read_csv(tsv_path, sep="\t")
    required_columns = ['Old ID', 'New ID']
    if not all(col in df.columns for col in required_columns):
        sys.exit(f"Error: TSV file must contain columns: {', '.join(required_columns)}")

    for _, row in df.iterrows():
        old_id = row['Old ID']
        new_id = row['New ID']

        try:
            org = Metadata.get(Metadata.short_name == old_id)
            if old_id == new_id:
                print(f"Skipping {old_id}: New ID is the same.")
                continue

            if dry_run:
                print(f"[Dry Run] Updating {old_id} -> {new_id}")
            else:
                # Update Metadata.short_name
                org.short_name = new_id
                org.save()

                # Update Sequences headers (name field)
                update = (
                    Sequences
                    .update({Sequences.name: fn.REPLACE(Sequences.name, old_id, new_id)})
                    .where(Sequences.organism == org)
                )
                updated_rows = update.execute()

                print(f"Updated {old_id} -> {new_id}")

        except Metadata.DoesNotExist:
            print(f"Skipping: {old_id} not found in metadata table.")
    
    if not dry_run:
        build_database.main(threads, no_og_file=True, threshold=0.1)


if __name__ == '__main__':
    description = 'Explore database'
    parser, optional, required = help_formatter.initialize_argparse(name='explore_database.py',
                                                                    desc=description,
                                                                    usage='explore_database.py [OPTIONS]')

    # Optional Arguments
    optional.add_argument('-d', '--database', type=str,
                          help=textwrap.dedent("""\
                          Path to the database folder if config.ini is not used.
                              """))
    optional.add_argument('-t', '--higher_taxonomy', action='store_true',
                          help="Show metadata grouped by higher taxonomy.")
    optional.add_argument('-l', '--lower_taxonomy', action='store_true',
                          help="Show metadata grouped by lower taxonomy.")
    optional.add_argument('-r', '--get_higher', default=None,
                          help="Return ortholog/paralog counts for a given Higher Taxonomy.")
    optional.add_argument('-w', '--get_lower', default=None,
                          help="Return ortholog/paralog counts for a given Lower Taxonomy.")
    optional.add_argument('-o', '--get_org', default=None,
                          help="Return information for a given organism by short name.")
    optional.add_argument('--update_metadata', type=str,
                          help="Update existing metadata rows using this TSV file. See documentation for format.")
    optional.add_argument('--show_taxonomy_colors', action='store_true',
                          help="Display all taxonomy terms with their associated colors.")
    optional.add_argument('--update_taxonomy_colors', type=str,
                          help="TSV file to update taxonomy colors. Columns: Taxonomy, Color.")
    optional.add_argument('--update_unique_ids', type=str,
                          help="TSV file to update Unique IDs and sequence headers. Columns: Old ID, New ID.")
    optional.add_argument('--threads', type=int,
                          help="Number of threads. Default is 1. Only to be used with --update_unique_ids.")
    optional.add_argument('--dry_run', default=False, action='store_true',
                          help="Do not update the database, just print what would be changed.\n" \
                          "For use with --update_metadata, --update_taxonomy_colors, or --update_unique_ids.")

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, out_dir=False, inp_dir=False)

    config = configparser.ConfigParser()
    config.read('config.ini')

    if args.database:
        dfo = str(Path(args.database).resolve())
    else:
        dfo = str(Path(config['PATHS']['database_folder']).resolve())

    # Connect to SQLite database using Peewee
    db_path = os.path.join(dfo, 'phylofisher.db')
    database.init(db_path) 
    database.connect()

    metadata_handle = MetadataHandler()

    if args.higher_taxonomy:
        print(metadata_handle.higher_taxonomy().to_string(index=False))

    elif args.lower_taxonomy:
        print(metadata_handle.lower_taxonomy().to_string(index=False))

    elif args.get_higher:
        print(metadata_handle.get_higher(args.get_higher).to_string(index=False))

    elif args.get_lower:
        print(metadata_handle.get_lower(args.get_lower).to_string(index=False))

    elif args.get_org:
        print(metadata_handle.get_org(args.get_org).to_string())
    
    elif args.show_taxonomy_colors:
        show_taxonomy_colors()
    
    elif args.update_taxonomy_colors:
        update_taxonomy_colors(args.update_taxonomy_colors, args.dry_run)

    elif args.update_metadata:
        update_metadata_table(args.update_metadata, args.dry_run)
    
    elif args.update_unique_ids:
        update_unique_ids(args.threads, args.update_unique_ids, args.dry_run)


    database.close()
