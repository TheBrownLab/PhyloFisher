#!/usr/bin/env python
import configparser
import os
import shutil
import sys
import tarfile
import textwrap
from datetime import date
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofisher import help_formatter, tools


def parse_collapse_tsv():
    '''
    Parses to collapse tsv file

    :return: chimera info dictionary with higher, lower, and unique ids
    :rtype: dict
    '''
    with open(args.input, 'r') as infile:
        collapse_dict = dict()
        for line in infile:
            # Check to make sure the correct number of columns are present
            try:
                chimera_id, long_name, higher, lower, to_collapse = line.strip().split('\t')
            except ValueError:
                raise Exception('Error parsing to_collapse.tsv. Please make sure the file is formatted correctly.')
            
            # Add parsed information to dictionary
            collapse_dict[chimera_id] = {
                'long_name': long_name,
                'higher': higher,
                'lower': lower,
                'to_collapse': to_collapse.split(',')
            }
    
    return collapse_dict


def check_metadata(collapse_dict):
    '''
    Checks to make sure taxa provided in to_collapse.tsv are actually in the metadata

    :param collapse_dict: information about taxa to collapse
    :type collapse_dict: dict
    '''
    to_collapse = [taxon for k in collapse_dict.keys() for taxon in collapse_dict[k]['to_collapse']]
    taxon_dict = dict.fromkeys(to_collapse, False)

    # Loop through metadata and check if taxa are in metadata
    with open(metadata, 'r') as infile:
        for line in infile:
            line = line.strip()
            short_name = line.split('\t')[0]
            if short_name in to_collapse:
                taxon_dict[short_name] = True

    # Check if all taxa are in metadata
    for taxon in taxon_dict.keys():
        if taxon_dict[taxon] is False:
            print(f'{taxon} is not in the dataset.')

    # Raise error if any taxa are not in metadata
    if False in taxon_dict.values():
        raise Exception('Taxa not collapsed. Please fix errors above.')


def add_to_meta(collapse_dict):
    '''
    Transfer input metadata for a given organism to dataset metadata.

    :param collapse_dict: information about taxa to collapse
    :type collapse_dict: dict
    '''
    for chimera_id in collapse_dict.keys():
        # Create new row dictionary
        new_row = {'Unique ID': chimera_id,
                'Long Name': collapse_dict[chimera_id]['long_name'],
                'Higher Taxonomy': collapse_dict[chimera_id]['higher'],
                'Lower Taxonomy': collapse_dict[chimera_id]['lower'],
                'Source': 'taxon_collapser.py',
                'Data Type': 'Collapsed'
                }

        # Read metadata into dataframe
        df = pd.read_csv(metadata, delimiter='\t')
        df.index = df['Unique ID']

        # Add new row to dataframe
        df = df.append(new_row, ignore_index=True)
        df = df.sort_values(by=['Higher Taxonomy', 'Lower Taxonomy', 'Unique ID'])
        df.to_csv(str(Path(dfo, 'metadata.tsv')), sep='\t', index=False)


def collapse(collapse_dict):
    '''
    Collapses sequences in the dataset.

    :param collapse_dict: information about taxa to collapse
    :type collapse_dict: dict
    '''

    # Collapse Ortholgogs
    ortho_dir = f'{dfo}/orthologs'
    gene_files = [os.path.join(ortho_dir, gene) for gene in os.listdir(ortho_dir) if gene.endswith('.fas')]
    for gene_file in gene_files:
        record_dict = {k: SeqRecord(Seq(''), id='', name='', description='') for k in collapse_dict.keys()}
        records = []
        with open(gene_file, 'r') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                for chimera_id in collapse_dict.keys():
                    if record.description in collapse_dict[chimera_id]['to_collapse']:
                        if len(record.seq) > len(record_dict[chimera_id].seq):
                            record_dict[chimera_id].id = chimera_id
                            record_dict[chimera_id].seq = record.seq

                # Append original sequences to records
                records.append(record)

        # Append collapsed sequences to records
        for record in record_dict.values():
            if str(record.seq) != '':
                records.append(record)

        # Write ortholog file post collapse
        with open(gene_file, 'w') as outfile:
            SeqIO.write(records, outfile, 'fasta')

        record_dict.clear()

    # Collapse Paralogs
    para_dir = f'{dfo}/paralogs'
    gene_files = [os.path.join(para_dir, gene) for gene in os.listdir(para_dir) if gene.endswith('.fas')]
    for gene_file in gene_files:
        records = []
        with open(gene_file, 'r') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                for chimera_id in collapse_dict.keys():
                    if record.description.split('..')[0] in collapse_dict[chimera_id]['to_collapse']:
                        records.append(SeqRecord(record.seq, id=f'{chimera_id}..{record.description.split("..")[1]}', name='', description=''))
                
                records.append(record)

        # Write ortholog file post collapse
        with open(gene_file, 'w') as outfile:
            SeqIO.write(records, outfile, 'fasta')

        record_dict.clear()

    add_to_meta(collapse_dict)


if __name__ == '__main__':
    description = 'Collapses dataset entries into a single taxon'
    parser, optional, required = help_formatter.initialize_argparse(name='taxon_collapser.py',
                                                                    desc=description,
                                                                    usage='taxon_collapser.py [OPTIONS] -i <taxa>')

    required.add_argument('-i', '--input', type=str, metavar='to_collapse.tsv', default=None,
                          help=textwrap.dedent("""\
                          A .tsv containing a Unique ID, long name, higher taxonomy, lower taxonomy, and a list of Unique IDs to collapse
                          for each chimera with one on each line.
                          
                          Example:
                            Chimera_ID 1\tLong Name\tHigher_Tax\tLower_Tax\tUnique_ID_1,Unique_ID_2
                            Chimera_ID 2\tLong Name\tHigher_Tax\tLower_Tax\tUnique_ID_3,Unique_ID_4
                        """))
    args = help_formatter.get_args(parser, optional, required,
                                   out_dir=False, pre_suf=False, inp_dir=False)

    # Config Parser
    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    metadata = str(os.path.join(dfo, 'metadata.tsv'))

    tools.backup(dfo)
    collapse_dict = parse_collapse_tsv()
    check_metadata(collapse_dict)
    collapse(collapse_dict)
