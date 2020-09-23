#!/usr/bin/env python
import configparser
import csv
import os
import textwrap
from glob import glob
from pathlib import Path

from Bio import SeqIO

from phylofisher import help_formatter


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


def parse_input(infile):
    file_contents = set()
    with open(infile, 'r') as infile:
        for line in infile:
            line = line.strip()
            file_contents.add(line)

    return list(file_contents)


def delete_group_org(orgs_=None, groups=None):
    meta = os.path.join(dfo, 'metadata.tsv')
    if groups:
        group_set = set(groups)
    else:
        group_set = set()

    if orgs_:
        orgs = set(orgs_)
    else:
        orgs = set()

    lines = parse_metadata()
    orgs_to_del = set()
    with open(meta, 'w') as out_file:
        res = csv.writer(out_file, delimiter='\t')
        for line in lines:
            if line[2] in group_set:
                orgs_to_del.add(line[0])
            elif line[3] in group_set:
                orgs_to_del.add(line[0])
            elif line[0] in orgs:
                orgs_to_del.add(line[0])
            else:
                res.writerow(line)
    delete_homologs(orgs_to_del)


# TODO input as a file
if __name__ == '__main__':
    description = 'Deletes taxa and/or taxonomic groups from the database'
    parser, optional, required = help_formatter.initialize_argparse(name='purge.py',
                                                                    desc=description,
                                                                    usage='purge.py [OPTIONS] -i <in_dir>')

    # Optional Arguments
    optional.add_argument('-o', '--orgs', type=str, metavar='groups.txt',
                          help=textwrap.dedent("""\
                              Path to text file containing Unique IDs of organisms for deletion.
                              Example:
                              UniqueID-1
                              UniqueID-2
                              """))
    optional.add_argument('-g', '--tax_groups', metavar='groups.txt', type=str,
                          help=textwrap.dedent("""\
                              Path to text file containing taxonomic groups for deletion.
                              Example:
                              Group1
                              Group2
                              """))

    in_help = 'Path to database directory'
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, in_help=in_help, out_dir=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    multi_input = os.path.abspath(config['PATHS']['input_file'])

    if args.orgs:
        orgs = parse_input(args.orgs)
        delete_group_org(orgs_=orgs)
    if args.tax_groups:
        tax_groups = parse_input(args.tax_groups)
        delete_group_org(groups=tax_groups)


