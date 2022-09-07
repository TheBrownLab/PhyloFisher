#!/usr/bin/env python
import configparser
import csv
import os
import shutil
import sys
import textwrap
from glob import glob
from pathlib import Path

from Bio import SeqIO

from phylofisher import help_formatter, tools


def parse_metadata():
    meta = os.path.join(dfo, 'metadata.tsv')
    with open(meta, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        lines = list(reader)
        return lines


def parse_input():
    file_contents = set()
    with open(args.input, 'r') as infile:
        for line in infile:
            line = line.strip()
            file_contents.add(line)

    return list(file_contents)


def check_metadata():
    collapsed_taxa = []
    for item in parse_input():
        for line in parse_metadata():
            if item not in line:
                sys.exit(f'{item} is not in the database. Please check your input file.')
            if 'taxon_collapser.py' in line:
                collapsed_taxa.append(line[0])

    return collapsed_taxa


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


def delete_proteomes(org_set, collapsed_taxa):
    for org in org_set:
        if org not in collapsed_taxa:
            os.remove(os.path.join(dfo, 'proteomes', f'{org}.faa.tar.gz'))


def purge(collapsed_taxa):
    to_remove = parse_input()
    meta = os.path.join(dfo, 'metadata.tsv')

    lines = parse_metadata()
    orgs_to_del = set()
    
    with open(meta, 'w') as out_file:
        res = csv.writer(out_file, delimiter='\t')
        for line in lines:
            if line[2] in to_remove:
                orgs_to_del.add(line[0])
            elif line[3] in to_remove:
                orgs_to_del.add(line[0])
            elif line[0] in to_remove:
                orgs_to_del.add(line[0])
            else:
                res.writerow(line)
                
    delete_homologs(orgs_to_del)
    delete_proteomes(orgs_to_del, collapsed_taxa)


# TODO input as a file
if __name__ == '__main__':
    description = 'Deletes taxa and/or taxonomic groups from the database'
    parser, optional, required = help_formatter.initialize_argparse(name='purge.py',
                                                                    desc=description,
                                                                    usage='purge.py [OPTIONS] -i to_purge.txt -d path/to/database')

    # Optional Arguments
    required.add_argument('-i', '--input', type=str, metavar='to_purge.txt',
                          help=textwrap.dedent("""\
                          Path to text file containing Unique IDs and Taxonomic designations of organisms for deletion.
                           """))
    required.add_argument('-d', '--database', metavar='<input_dir>', type=str,
                          help=textwrap.dedent("""\
                          Path to database to purge.
                          """))

    in_help = 'Path to database directory'
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    dfo = os.path.abspath(args.database)

    collapsed_taxa = check_metadata()
    tools.backup(dfo)
    purge(collapsed_taxa)
