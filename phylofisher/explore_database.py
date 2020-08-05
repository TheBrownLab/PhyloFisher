#!/usr/bin/env python
import os
import configparser
from glob import glob
from pathlib import Path
from Bio import SeqIO
import pandas as pd
from phylofisher import help_formatter
import textwrap

pd.options.display.float_format = '{:,.0f}'.format


class Metadata:
    """
    Class which store and combine information from metadata.tsv and 
    orthologs and paralogs from database.
    """
    def __init__(self, metadata, ortholog_folder, paralog_folder):
        self.metadata = pd.read_csv(metadata, sep='\t')
        self.metadata = self.metadata[['Unique ID', 'Long Name', 
                                    'Higher Taxonomy', 'Lower Taxonomy']]
        self.ortholog_folder = ortholog_folder
        self.paralog_folder = paralog_folder
        self.metadata['orthologs'] =  self.metadata['Unique ID'].map(self.parse_orthologs())
        self.metadata['paralogs'] =  self.metadata['Unique ID'].map(self.parse_paralogs())

    def parse_orthologs(self):
        """
        Count orthologs for all organisms.
        return: dictionary[org] = count
        """
        counted_orthologs = {}
        for file in glob(f'{self.ortholog_folder}/*.fas'):
            for record in SeqIO.parse(file, 'fasta'):
                if record.name not in counted_orthologs:
                    counted_orthologs[record.name] = 1
                else:
                    counted_orthologs[record.name] += 1
        return counted_orthologs

    def parse_paralogs(self):
        """
        Count paralogs for all organisms.
        return: dictionary[org] = count
        """
        counted_paralogs = {}
        for file in glob(f'{self.paralog_folder}/*.fas'):
            for record in SeqIO.parse(file, 'fasta'):
                name = record.name.split('.')[0]
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
        ).size().reset_index().rename(columns={0:'organisms'})


    def lower_taxonomy(self):
        """
        Return table with metadata grouped according to a lower taxonomy.
        """
        return self.metadata.groupby(['Higher Taxonomy','Lower Taxonomy']
        ).size().reset_index().rename(columns={0:'organisms'})


    def get_higher(self, term):
        """
        Return table with orthologs, paralogs for organisms from a given Higher Taxonomy.
        """
        return self.metadata[self.metadata['Higher Taxonomy'
        ] == term].sort_values(['Lower Taxonomy', 'Unique ID']).reset_index(drop=True)


    def get_lower(self, term):
        """
        Return table with orthologs, paralogs for organisms from a given Lower Taxonomy.
        """
        return self.metadata[self.metadata['Lower Taxonomy'
        ] == term].sort_values('Unique ID').reset_index(drop=True)




if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    ortholog_folder  = str(Path(dfo, f'orthologs'))
    paralog_folder = str(Path(dfo, f'paralogs'))
    metadata = str(os.path.join(dfo, 'metadata.tsv'))
    metadata_handle = Metadata(metadata, ortholog_folder, paralog_folder)
    

    description = 'Explore database'
    parser, optional, required = help_formatter.initialize_argparse(name='explore_database.py',
                                                                    desc=description,
                                                                    usage='explore_database.py [OPTIONS]')

    # Optional Arguments
    optional.add_argument('-t','--higher_taxonomy', action='store_true',
                          help=textwrap.dedent("""\
                          Show table with metadata grouped according to the higher taxonomy.
                              """))

    optional.add_argument('-l','--lower_taxonomy', action='store_true',
                          help=textwrap.dedent("""\
                          Show table with metadata grouped according to the lower taxonomy.
                              """))

    optional.add_argument('-r','--get_higher', default=None,
                          help=textwrap.dedent("""\
                          Return table with orthologs, paralogs for organisms from a given Higher Taxonomy.
                              """))

    optional.add_argument('-w','--get_lower', default=None,
                          help=textwrap.dedent("""\
                          Return table with orthologs, paralogs for organisms from a given Lower Taxonomy.
                              """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, out_dir=False, inp_dir=False)

    if args.higher_taxonomy:
        print(metadata_handle.higher_taxonomy().to_string())

    elif args.lower_taxonomy:
        print(metadata_handle.lower_taxonomy().to_string())

    elif args.get_higher:
        print(metadata_handle.get_higher(args.get_higher).to_string())

    elif args.get_lower:
        print(metadata_handle.get_lower(args.get_lower).to_string())