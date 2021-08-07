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
    with open(args.input, 'r') as infile:
        collapse_dict = dict()
        for line in infile:
            collapse_dict[line.strip().split('\t')[0]] = line.strip().split('\t')[1:]
    return collapse_dict


def check_metadata():
    """
    Checks to make sure taxa provided in to_collapse.csv are actually in the meta data
    """
    with open(args.input, 'r') as infile:
        to_collapse = set()
        for taxon in [line.strip().split(',')[1:] for line in infile]:
            to_collapse.update(taxon)

    taxon_dict = dict.fromkeys(to_collapse, False)

    with open(metadata, 'r') as infile:
        for line in infile:
            line = line.strip()

            short_name = line.split('\t')[0]
            if short_name in to_collapse:
                taxon_dict[short_name] = True

    for taxon in taxon_dict.keys():
        if taxon_dict[taxon] is False:
            print(f'{taxon} is not in the dataset.')

    if False in taxon_dict.values():
        sys.exit('Taxa not collapsed. Please fix errors above.')


def add_to_meta(abbrev, input_info):
    """"Transfer input metadata for a given organism
    to dataset metadata.
    input:  short name of organism
    return: None
    """

    new_row = {'Unique ID': abbrev,
               'Long Name': input_info[abbrev][0],
               'Higher Taxonomy': input_info[abbrev][1],
               'Lower Taxonomy': input_info[abbrev][2],
               'Source': 'taxon_collapser.py',
               'Data Type': 'Collapsed'
               }

    df = pd.read_csv(metadata, delimiter='\t')
    df.index = df['Unique ID']

    df = df.append(new_row, ignore_index=True)
    df = df.sort_values(by=['Higher Taxonomy', 'Lower Taxonomy', 'Unique ID'])
    df.to_csv(str(Path(dfo, 'metadata.tsv')), sep='\t', index=False)


def collapse():
    collapse_dict = parse_collapse_tsv()

    # Collapse Ortholgogs
    ortho_dir = f'{dfo}/orthologs'
    gene_files = [os.path.join(ortho_dir, gene) for gene in os.listdir(ortho_dir) if gene.endswith('.fas')]
    for gene_file in gene_files:
        record_dict = {k: SeqRecord(Seq(''), id='', name='', description='') for k in collapse_dict.keys()}
        records = []
        with open(gene_file, 'r') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                for key in collapse_dict.keys():
                    if record.description in collapse_dict[key]:
                        if len(record.seq) > len(record_dict[key].seq):
                            record_dict[key].id = key
                            record_dict[key].seq = record.seq

                records.append(record)

        for record in record_dict.values():
            if str(record.seq) != '':
                records.append(record)

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
                for key in collapse_dict.keys():
                    if record.description.split('..')[0] in collapse_dict[key]:
                        records.append(SeqRecord(record.seq, id=f'{key}..{record.description.split("..")[1]}', name='', description=''))
                
                records.append(record)

        with open(gene_file, 'w') as outfile:
            SeqIO.write(records, outfile, 'fasta')

        record_dict.clear()

    for key in collapse_dict.keys():
        add_to_meta(key, collapse_dict)


if __name__ == '__main__':
    description = 'Collapses dataset entries into a single taxon'
    parser, optional, required = help_formatter.initialize_argparse(name='taxon_collapser.py',
                                                                    desc=description,
                                                                    usage='taxon_collapser.py [OPTIONS] -i <taxa>')

    required.add_argument('-i', '--input', type=str, metavar='to_collapse.tsv', default=None,
                          help=textwrap.dedent("""\
                          A .tsv containing a Unique ID, higher and lower taxonomic designations, 
                           and the Unique IDs of the taxa to collapse, for each chimera one per line"""))
    args = help_formatter.get_args(parser, optional, required,
                                   out_dir=False, pre_suf=False, inp_dir=False)

    # Config Parser
    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    metadata = str(os.path.join(dfo, 'metadata.tsv'))

    tools.backup(dfo)
    check_metadata()
    collapse()
