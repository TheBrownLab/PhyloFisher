#!/usr/bin/env python
import configparser
import os, sys
import textwrap
from glob import glob
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from phylofisher import help_formatter

pd.options.display.float_format = '{:,.0f}'.format


class Metadata:
    """
    Class which store and combine information from metadata.tsv and 
    orthologs and paralogs from database.
    """

    def __init__(self, metadata, ortholog_folder, paralog_folder):
        self.metadata = pd.read_csv(metadata, sep='\t')

        self.ortholog_folder = ortholog_folder
        self.paralog_folder = paralog_folder
        self.metadata['Orthologs'] = self.metadata['Unique ID'].map(self.parse_orthologs())
        self.metadata['Paralogs'] = self.metadata['Unique ID'].map(self.parse_paralogs())

        self.org_meta = self.metadata[['Unique ID', 'Long Name', 'Higher Taxonomy', 'Lower Taxonomy',
                                       'Data Type', 'Orthologs', 'Paralogs', 'Source']]
        self.metadata = self.metadata[['Unique ID', 'Long Name', 'Higher Taxonomy', 'Lower Taxonomy',
                                       'Data Type', 'Orthologs', 'Paralogs']]


    def parse_orthologs(self):
        """
        Count orthologs for all organisms.
        return: dictionary[org] = count
        """
        counted_orthologs = {}
        for file in glob(f'{self.ortholog_folder}/*.fas'):
            for record in SeqIO.parse(file, 'fasta'):
                if record.description not in counted_orthologs:
                    counted_orthologs[record.description] = 1
                else:
                    counted_orthologs[record.description] += 1
        return counted_orthologs

    def parse_paralogs(self):
        """
        Count paralogs for all organisms.
        return: dictionary[org] = count
        """
        counted_paralogs = {}
        for file in glob(f'{self.paralog_folder}/*.fas'):
            for record in SeqIO.parse(file, 'fasta'):
                name = record.description.split('.')[0]
                if name not in counted_paralogs:
                    counted_paralogs[name] = 1
                else:
                    counted_paralogs[name] += 1
        return counted_paralogs

    def higher_taxonomy(self):
        """
        Return table with metadata grouped according to a higher taxonomy.
        """
        return self.metadata.groupby(['Higher Taxonomy']
                                     ).size().reset_index().rename(columns={0: 'Organisms'})

    def lower_taxonomy(self):
        """
        Return table with metadata grouped according to a lower taxonomy.
        """
        return self.metadata.groupby(['Higher Taxonomy', 'Lower Taxonomy']
                                     ).size().reset_index().rename(columns={0: 'Organisms'})

    def get_higher(self, term):
        """
        Return table with orthologs, paralogs for organisms from a given Higher Taxonomy.
        """
        try:
            ret = self.metadata[self.metadata['Higher Taxonomy'] == term].sort_values(['Lower Taxonomy', 'Unique ID']).reset_index(drop=True)
        except KeyError:
            sys.exit(f'{term} is not a valid Higher Taxonomy')
        return ret

    def get_lower(self, term):
        """
        Return table with orthologs, paralogs for organisms from a given Lower Taxonomy.
        """
        try:
            ret = self.metadata[self.metadata['Lower Taxonomy'] == term].sort_values('Unique ID').reset_index(drop=True)
        except KeyError:
            sys.exit(f'{term} is not a valid Lower Taxonomy')
        return ret 

    def get_org(self, term):
        """
        Return ifnormation about a given organism using its short name.
        """
        org_df = self.org_meta.copy()
        org_df = org_df.set_index('Unique ID')
        try:
            df = pd.DataFrame([org_df.loc[term, :]]).transpose()
        except KeyError:
            sys.exit(f'{term} is not a valid UniqueID')
        return df


if __name__ == '__main__':
    description = 'Explore database'
    parser, optional, required = help_formatter.initialize_argparse(name='explore_database.py',
                                                                    desc=description,
                                                                    usage='explore_database.py [OPTIONS]')

    # Optional Arguments
    optional.add_argument('-d', '--database', type=str,  # metvar='<path/to/database>',
                          help=textwrap.dedent("""\
                          Path to the database folder if config.ini is not used.
                              """))
    optional.add_argument('-t', '--higher_taxonomy', action='store_true',
                          help=textwrap.dedent("""\
                          Show table with metadata grouped according to the higher taxonomy.
                              """))

    optional.add_argument('-l', '--lower_taxonomy', action='store_true',
                          help=textwrap.dedent("""\
                          Show table with metadata grouped according to the lower taxonomy.
                              """))

    optional.add_argument('-r', '--get_higher', default=None,
                          help=textwrap.dedent("""\
                          Return table with orthologs, paralogs for organisms from a given Higher Taxonomy.
                              """))

    optional.add_argument('-w', '--get_lower', default=None,
                          help=textwrap.dedent("""\
                          Return table with orthologs, paralogs for organisms from a given Lower Taxonomy.
                              """))

    optional.add_argument('-o', '--get_org', default=None,
                          help=textwrap.dedent("""\
                          Return information about a given organism using its short name.
                              """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, out_dir=False, inp_dir=False)

    config = configparser.ConfigParser()
    config.read('config.ini')

    if args.database:
        dfo = str(Path(args.database).resolve())
    else:
        dfo = str(Path(config['PATHS']['database_folder']).resolve())

    ortholog_folder = str(Path(dfo, f'orthologs'))
    paralog_folder = str(Path(dfo, f'paralogs'))
    metadata = str(os.path.join(dfo, 'metadata.tsv'))
    metadata_handle = Metadata(metadata, ortholog_folder, paralog_folder)

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
