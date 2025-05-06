#!/usr/bin/env python
import configparser
from email.policy import default
import glob
import os
import random
import shutil
import string
import sys
import tarfile
import textwrap
from collections import defaultdict
from datetime import date
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from phylofisher import help_formatter, tools
from phylofisher.utilities import build_database
from phylofisher.db_map import database, Genes, Taxonomies, Metadata, Sequences


class UnknownStatusError(Exception):
    '''
    Exception Raised for unknown status
    '''
    def __init__(self, unique_id, table):
        self.unique_id = unique_id
        self.table = table
        self.message =  f'Unknown status for {self.unique_id} in {self.table}'
        super().__init__(self.message)


def taxa_to_exclude():
    '''
    Parse taxa to exclude from dataset addition.

    :return: list of unique IDs to exclude
    :rtype: list
    '''
    to_skip = []
    with open(args.to_exclude, 'r') as infile:
        for line in infile:
            unique_id = line.strip().split(',')[0]
            to_skip.append(unique_id)

    return to_skip


def parse_input(input_metadata, orgs_to_exc):
    '''
    Parse input metadata.

    :param input_metadata: path to input metadata file
    :type input_metadata: str
    :param orgs_to_exc: list of unique IDs to exclude from dataset addition
    :type orgs_to_exc: list
    :return: dictionary of input metadata
    :rtype: dict
    '''
    orgs_to_exc = set()

    input_info = defaultdict(dict)
    with open(input_metadata, 'r') as infile:
        for line in infile:
            line = line.strip()
            sline = line.split('\t')
            abbrev = sline[2]

            if 'Location' in line:
                pass
            elif abbrev in orgs_to_exc:
                pass
            else:
                group = sline[3].strip()
                full_name = sline[6].strip()
                subtax = sline[4]
                data_type = sline[7]
                notes = sline[8].strip()
                input_info[abbrev]['tax'] = group
                input_info[abbrev]['full_name'] = full_name
                input_info[abbrev]['subtax'] = subtax
                input_info[abbrev]['data_type'] = data_type
                input_info[abbrev]['notes'] = notes

    return input_info


def collect_seqs(gene):
    '''
    Collect all sequences for a given gene from fasta folder (fisher result) and store them in a dictionary.

    :param gene: name of the gene
    :type gene: str
    :return: dictionary of all sequences for a given gene
    :rtype: dict
    '''
    seq_dict = {}
    for record in SeqIO.parse(f'{fisher_dir}/{gene}.fas', 'fasta'):
        if record.name.count('_') == 3:
            abbrev, _, _, quality = record.name.split('_')
            qname = f'{abbrev}_{quality}'
            seq_dict[qname] = record
        else:
            seq_dict[record.name] = record
    return seq_dict


def parse_table(table):
    '''
    Collect information what has to be changed in the PhyloFisher dataset (for a given gene) 
    according to information from parsed single gene tree.

    :param table: path to parsed tsv file of one gene
    :type table: str
    :return: dictionaries of new orthologs and paralogs
    :rtype: dict
    '''

    # gene name
    gene = table.split('/')[-1].split('_parsed')[0]
    g_id = Genes.select(Genes.id).where(Genes.name == gene).get().id

    # parse orthologs for a given gene into a dictionary
    orthologs = {}
    db_guery = Sequences.select(
        Sequences.header, Sequences.sequence).where(
            Sequences.is_paralog == False, Sequences.gene_id == g_id)
    for q in db_guery:
        orthologs[q.header] = q.sequence
    
    
    # parse paralogs into a dictionary (if they exist); else make an empty dict
    paralogs = {}
    db_guery = Sequences.select(
        Sequences.header, Sequences.sequence, Sequences.id).where(
            Sequences.is_paralog == True, Sequences.gene_id == g_id)
    for q in db_guery:
        paralogs[f'{q.header}..p{q.id}'] = q.sequence

    # collect all candidate sequences for a given gen from the fisher.py result
    seq_dict = collect_seqs(gene)

    for line in open(table):
        # for orthologs from the dataset
        tree_name, tax, status = line.split('\t')
        status = status.strip()  # o,p,d (ortholog, paralog, delete)
        abbrev = tree_name.split('@')[-1]
        if tree_name.count('_') != 3 and '..' not in abbrev:
            record = seq_dict[abbrev]
            db_guery = Sequences.select(
                Sequences.id, Sequences.header, Sequences.sequence).where(
                    Sequences.is_paralog == False, Sequences.gene_id == g_id, Sequences.header == abbrev)
            assert db_guery.count() == 1
            if status == 'd':
                # delete sequence from orthologs
                del orthologs[abbrev]
            elif status == 'p':
                # chance status from paralog to ortholog
                # prepare paralog name
                pname = f'{abbrev}..p{db_guery.get().id}'
                paralogs[pname] = record
                # delete sequence from orthologs
                del orthologs[abbrev]

    for line in open(table):
        # for new sequences
        tree_name, tax, status = line.split('\t')
        status = status.strip()
        abbrev = tree_name.split('@')[-1]
        if tree_name.count('_') == 3:
            quality = tree_name.split('_')[2]
            qname = f'{abbrev}_{quality}'
            record = seq_dict[qname]
            if status == 'o':
                # add to orthologs
                orthologs[abbrev] = record.seq
            elif status == 'p':
                # add to paralogs
                pname = f'{abbrev}..p{db_guery.get().id}'
                paralogs[pname] = record.seq

    for line in open(table):
        # for paralogs from the dataset
        tree_name, tax, status = line.split('\t')
        status = status.strip()
        name = tree_name.split('@')[-1]
        abbrev = name.split('.')[0]
        db_guery = Sequences.select(
            Sequences.id, Sequences.header, Sequences.sequence).where(
                Sequences.is_paralog == True, Sequences.gene_id == g_id, Sequences.header == abbrev)
        if '..' in name:
            if status == 'o':
                record = paralogs[name]
                # for cases when ortholog has not survived trimming
                if abbrev in orthologs:
                    # prepare paralog name
                    pname = f'{abbrev}..p{db_guery.get().id}'
                    # chance status from paralog to ortholog
                    paralogs[pname] = orthologs[abbrev]
                    # record in orthologs will be replaced
                    # with a new one in the next step
                # add sequence to orthologs 
                orthologs[abbrev] = record
                # delete sequence from paralogs
                del paralogs[name]

            
    return orthologs, paralogs


def add_to_meta(abbrev):
    '''
    Transfer input metadata for a given organism to database metadata.

    :param abbrev: abbreviation of the organism
    :type abbrev: str
    '''
    tax = input_info[abbrev]['tax']
    subtax = input_info[abbrev]['subtax']
    
    # Check if taxonomies exist in the database
    if '*' in tax:
        tax = tax.replace('*', '')
        Taxonomies.insert(taxonomy=tax).execute()
    if '*' in subtax:
        subtax = subtax.replace('*', '')
        Taxonomies.insert(taxonomy=subtax).execute()
    
    # Insert new metadata rows
    q = Metadata.insert(short_name=abbrev,
                        long_name=input_info[abbrev]['full_name'],
                        higher_taxonomy_id=Taxonomies.select().where(Taxonomies.taxonomy == tax).get().id,
                        lower_taxonomy_id=Taxonomies.select().where(Taxonomies.taxonomy == subtax).get().id,
                        data_type=input_info[abbrev]['data_type'],
                        source=input_info[abbrev]['notes'])
    q.execute()


def new_database(table):
    '''
    Make changes in the PhyloFisher database according to information from parsed single homolog tree.

    :param table: path to parsed tsv file of one gene
    :type table: str
    '''
    gene = table.split('/')[-1].split('_parsed')[0]
    # Orthologs for a given gene
    orthologs_path = str(Path(dfo, f'orthologs/{gene}.fas'))
    # Paralogs for a given gene
    paralogs_path = str(Path(dfo, f'paralogs/{gene}_paralogs.fas'))

    orthologs, paralogs = parse_table(table)

    # make changes in orthologs
    for name, seq in orthologs.items():
        if name not in meta_orgs and name not in to_exclude:
            meta_orgs.add(name)
            add_to_meta(name)
        if name.split('..')[0] not in to_exclude:
            q = Sequences.insert(header=name,
                                 sequence=seq,
                                 is_paralog=False,
                                 gene_id=Genes.select().where(Genes.name == gene).get().id,
                                 metadata_id=Metadata.select().where(Metadata.short_name == name).get().id)
            q.execute()

    # make changes in paralogs
    for name, seq in paralogs.items():
        if name.split('..')[0] not in meta_orgs and name.split('..')[0] not in to_exclude:
            meta_orgs.add(name.split('.')[0])
            add_to_meta(name.split('.')[0])
        if name.split('..')[0] not in to_exclude:
            name = name.split('..')[0]
            q = Sequences.insert(header=name,
                                 sequence=seq,
                                 is_paralog=True,
                                 gene_id=Genes.select().where(Genes.name == gene).get().id,
                                 metadata_id=Metadata.select().where(Metadata.short_name == name).get().id)
            q.execute()


def rebuild_db():
    '''
    Rebuild the database after adding new data.
    '''
    cwd = os.getcwd()
    os.chdir(dfo)
    args.rename = None
    build_database.main(args, args.threads, True, 0.1)
    os.chdir(cwd)


def cp_proteomes():
    '''
    Copy proteomes from input metadata to database folder.
    '''
    with open(input_metadata, 'r') as infile:
        infile.readline()
        file_dict = {}
        for line in infile:
            s_line = line.strip().split('\t')
            file_dict[s_line[2]] = os.path.abspath(os.path.join(s_line[0], s_line[1]))

    # If tmp dir exists, remove it
    if os.path.isdir('tmp'):
        shutil.rmtree('tmp')
    
    os.mkdir('tmp')
    os.chdir('tmp')
    for key in file_dict.keys():
        # Check to make sure input proteome exits
        if not os.path.isfile(file_dict[key]):
            sys.exit(f'File {file_dict[key]} does not exist')
        with tarfile.open(f'{key}.faa.tar.gz', "x:gz") as tar:
            tar.add(file_dict[key], arcname=f'{key}.faa')
        shutil.copy(f'{key}.faa.tar.gz', f'{dfo}/proteomes/{key}.faa.tar.gz')

    os.chdir('..')
    shutil.rmtree('tmp')


def main():
    '''
    Loop through all parsed tables and add new data to the database.
    '''
    for table in glob.glob(f'{args.input}/*_parsed.tsv'):
        new_database(table)

    rebuild_db()
    cp_proteomes()


if __name__ == '__main__':
    desc = ('Apply parsing decisions and add new data to the database. \n\n'
            'NOTE: If apply_to_db.py fails for any reason, see backup_resoration.py \n'
            'and restore you database from the proper backup')
    parser, optional, required = help_formatter.initialize_argparse(name='apply_to_db.py',
                                                                    desc=desc,
                                                                    usage="apply_to_db.py -i <in_dir>")

    required.add_argument('-fi', '--fisher_dir', type=str, metavar='<fisher_dir>', required=True,
                          help=textwrap.dedent("""\
                          Path to fisher output directory to use for dataset addition.
                          """))

    optional.add_argument('--to_exclude', type=str, metavar='to_exclude.txt', default=None,
                          help=textwrap.dedent("""\
                          Path to .txt file containing Unique IDs of taxa to exclude from dataset 
                          addition with one taxon per line.
                          Example:
                            Unique ID (taxon 1)
                            Unique ID (taxon 2)
                            """))

    optional.add_argument('-t', '--threads', type=int, metavar='N', default=1,
                          help=textwrap.dedent("""\
                          Number of threads, where N is an integer.
                          Default: 1
                          """))

    in_help = 'Path to forest output directory to use for database addition.'
    args = help_formatter.get_args(parser, optional, required, out_dir=False, pre_suf=False, in_help=in_help)

    fisher_dir = args.fisher_dir

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    database.init(str(Path(dfo, 'phylofisher.db')))
    database.connect()

    tools.backup(dfo)
    
    if args.to_exclude:
        to_exclude = taxa_to_exclude()
    else:
        to_exclude = []
        
    input_metadata = os.path.abspath(config['PATHS']['input_file'])
    metadata = str(Path(dfo, 'metadata.tsv'))
    meta_orgs = set([x.short_name for x in Metadata.select(Metadata.short_name)])
    input_info = parse_input(input_metadata, to_exclude)
    main()

    database.close()