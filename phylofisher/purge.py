#!/usr/bin/env python
import os
import configparser
import argparse
import csv
import textwrap
import phylofisher.help_formatter
from Bio import SeqIO
from glob import glob
from pathlib import Path


def fasta_cleaner(file, org_set):
    records = list(SeqIO.parse(file, 'fasta'))
    with open(file, 'w') as res:
        for record in records:
            if record.name.split('.')[0] not in org_set:
                res.write(f'>{record.name}\n{record.seq}\n')


def delete_homologs(org_set):
    for folder in ['orthologs', 'paralogs']:
        files = glob(os.path.join(dfo, folder) + '/*.fas')
        for file in files:
            fasta_cleaner(file, org_set)


def parse_metadata():
    meta = os.path.join(dfo, 'metadata.tsv')
    with open(meta, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        lines = list(reader)
        return lines


def delete_group_org(orgs_=None, groups=None):
    meta = os.path.join(dfo, 'metadata.tsv')
    if groups:
        group_set = set(groups.split(','))
    else:
        group_set = set()

    if orgs_:
        orgs = set(orgs_.split(','))
    else:
        orgs = set()

    lines = parse_metadata()
    orgs_to_del = set()
    with open(meta, 'w') as out_file:
        res = csv.writer(out_file, delimiter='\t')
        for line in lines:
            if line[2] in group_set:
                orgs_to_del.add(line[0])
            elif line[0] in orgs:
                orgs_to_del.add(line[0])
            else:
                res.writerow(line)
    delete_homologs(orgs_to_del)


# TODO input as a file
if __name__ == '__main__':
    formatter = lambda prog: phylofisher.help_formatter.MyHelpFormatter(prog, max_help_position=100)
    description = 'Inspects single gene trees for contamination.'
    parser, optional, required = help_formatter.initialize_argparse(name='purge.py',
                                                                    desc=description,
                                                                    usage='forest.py [OPTIONS] -i <in_dir>')
    parser = argparse.ArgumentParser(prog='purge.py',
                                     description='Script for deleting orgs/taxonomic groups from the dataset',
                                     usage='purge.py [OPTIONS]',
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent("""\
                                         additional information:
                                            stuff
                                            """))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    # Required Arguments
    required.add_argument('-i', '--input', required=True, type=str, metavar='path/to/input/',
                          help=textwrap.dedent("""\
                              Path to input directory
                              """))

    # Optional Arguments
    optional.add_argument('-o', '--orgs', default="output", type=str, metavar='',
                          help=textwrap.dedent("""\
                              Short names of organisms for deletion: Org1,Org2,Org3
                              """))
    optional.add_argument('-g', '--tax_groups', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                              Names of taxonomic groups for deletion: Group1,Group2,Group3
                              """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                              Show this help message and exit.
                              """))

    parser._action_groups.append(optional)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    multi_input = os.path.abspath(config['PATHS']['input_file'])

    delete_group_org(args.orgs, args.tax_groups)

