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
from peewee import *

from phylofisher import help_formatter, tools
from phylofisher.db_map import database, Taxonomies, Metadata, Sequences


def parse_metadata():
    '''
    Queries metadata from SQLite database

    :return: list of tuples (short_name, higher_taxonomy, lower_taxonomy, long_name, source)
    :rtype: list
    '''
    higher = Taxonomies.alias('higher')
    lower = Taxonomies.alias('lower')
    
    query = (
        Metadata
        .select(
            Metadata.short_name,
            higher.taxonomy.alias('higher_taxonomy'),
            lower.taxonomy.alias('lower_taxonomy'),
            Metadata.long_name,
            Metadata.source
        )
        .join(higher, on=(Metadata.higher_taxonomy == higher.id))
        .switch(Metadata)
        .join(lower, on=(Metadata.lower_taxonomy == lower.id))
    )
    
    lines = []
    for row in query.dicts():
        lines.append([
            row['short_name'],
            row['long_name'],
            row['higher_taxonomy'],
            row['lower_taxonomy'],
            row['source']
        ])
    
    return lines


def parse_input():
    '''
    Parses input file containing taxa to remove

    :return: file contents
    :rtype: list
    '''
    file_contents = set()
    with open(args.input, 'r') as infile:
        for line in infile:
            line = line.strip()
            file_contents.add(line)

    return list(file_contents)


def check_metadata():
    '''
    Checks that taxa to remove are in database and returns a list of collapsed taxa

    :return: collapsed taxa
    :rtype: list
    '''
    collapsed_taxa = []
    all_metadata = []
    
    for line in parse_metadata():
        all_metadata += line
        if 'taxon_collapser.py' in line:
            collapsed_taxa.append(line[0])

    for item in parse_input():
        if item not in all_metadata:
            sys.exit(f'{item} is not in the database. Please check your input file.')

    return collapsed_taxa


def fasta_cleaner(file, org_set):
    '''
    Removes sequences from fasta file

    :param file: file to clean
    :type file: str
    :param org_set: organisms to remove
    :type org_set: set
    '''
    records = list(SeqIO.parse(file, 'fasta'))
    with open(file, 'w') as res:
        for record in records:
            if record.name.split('.')[0] not in org_set:
                res.write(f'>{record.name}\n{record.seq}\n')


def delete_homologs(org_set):
    '''
    Purges homologs from orthologs and paralogs directories

    :param org_set: organisms to remove
    :type org_set: set
    '''
    for folder in ['orthologs', 'paralogs']:
        files = glob(os.path.join(dfo, folder) + '/*.fas')
        for file in files:
            fasta_cleaner(file, org_set)


def purge(collapsed_taxa):
    '''
    Purges taxa from SQLite database

    :param collapsed_taxa: collapsed taxa
    :type collapsed_taxa: list
    '''
    to_remove = parse_input()
    lines = parse_metadata()
    orgs_to_del = set()
    metadata_to_del = []
    
    # Identify organisms to delete
    for line in lines:
        # line format: [short_name, long_name, higher_taxonomy, lower_taxonomy, source]
        if line[2] in to_remove:  # higher_taxonomy
            orgs_to_del.add(line[0])
            metadata_to_del.append(line[0])
        elif line[3] in to_remove:  # lower_taxonomy
            orgs_to_del.add(line[0])
            metadata_to_del.append(line[0])
        elif line[0] in to_remove:  # short_name
            orgs_to_del.add(line[0])
            metadata_to_del.append(line[0])
    
    # Delete sequences first (due to foreign key constraints)
    for org_name in metadata_to_del:
        meta = Metadata.get(Metadata.short_name == org_name)
        Sequences.delete().where(Sequences.metadata == meta).execute()
    
    # Delete metadata entries
    for org_name in metadata_to_del:
        Metadata.delete().where(Metadata.short_name == org_name).execute()
    
    # Also delete from fasta files
    delete_homologs(orgs_to_del)


if __name__ == '__main__':
    description = 'Deletes taxa and/or taxonomic groups from the database'
    parser, optional, required = help_formatter.initialize_argparse(name='purge.py',
                                                                    desc=description,
                                                                    usage='purge.py [OPTIONS] -i to_purge.txt -d path/to/database')

    # Required Arguments
    required.add_argument('-i', '--input', type=str, metavar='to_purge.txt',
                          help=textwrap.dedent("""\
                          Path to text file containing Unique IDs and Taxonomic designations of organisms for deletion.
                           """))
    required.add_argument('-d', '--database', metavar='<input_dir>', type=str,
                          help=textwrap.dedent("""\
                          Path to database to purge.
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    dfo = os.path.abspath(args.database)
    
    # Connect to SQLite database
    db_path = os.path.join(dfo, 'phylofisher.db')
    database.init(db_path)
    database.connect()

    collapsed_taxa = check_metadata()
    tools.backup(dfo)
    purge(collapsed_taxa)
    
    database.close()
