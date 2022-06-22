#!/usr/bin/env python
import configparser
from email.policy import default
import glob
import os
import random
import shutil
import string
import tarfile
import textwrap
from collections import defaultdict
from datetime import date
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from phylofisher import help_formatter, tools
from phylofisher.utilities import build_database


class UnknownStatusError(Exception):
    '''
    Exception Raised for unknown status
    '''
    def __init__(self, unique_id, table):
        self.unique_id = unique_id
        self.table = table
        self.message =  f'Unknown status for {self.unique_id} in {self.table}'
        super().__init__(self.message)


def dataset_orgs():
    """"Collect all short names from dataset metadata table."""
    orgs = set()
    with open(metadata) as md:
        # skip header
        next(md)
        for line_ in md:
            sline = line_.split('\t')
            abbrev = sline[0].strip()
            orgs.add(abbrev)
    return orgs


def taxa_to_exclude():
    to_skip = []
    with open(args.to_exclude, 'r') as infile:
        for line in infile:
            unique_id = line.strip().split(',')[0]
            to_skip.append(unique_id)

    return to_skip


def parse_input(input_metadata, orgs_to_exc):
    """"Parse input metadata.
    input: input metadata csv
    return: dictionary with info about input metadata
    """
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
    """Collect all sequences for a given gene from fasta folder (fisher result)
    and store them in a dictionary.
    input: gene name
    return: dictionary of all seqs for a given gene
    """
    seq_dict = {}
    for record in SeqIO.parse(f'{fisher_dir}/{gene}.fas', 'fasta'):
        if record.name.count('_') == 3:
            abbrev, _, _, quality = record.name.split('_')
            qname = f'{abbrev}_{quality}'
            seq_dict[qname] = record
        else:
            seq_dict[record.name] = record
    return seq_dict


def id_generator(size=5, chars=string.digits):
    """"Generate random number with 5 digits."""
    return ''.join(random.choice(chars) for _ in range(size))


def paralog_name(abbrev, keys):
    """"Prepare paralog name (short name + 5 random digits).
    Recursive function.
    example: Homosap..12345
    input: short name of an organism, names of already existing paralogs
    for a given organism
    return: unique paralog name"""
    id_ = id_generator()
    pname = f'{abbrev}..p{id_}'
    if pname not in keys:
        return pname
    else:
        paralog_name(abbrev, keys)


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

    # paths to orthologs and in the dataset folder
    orthologs_path = str(Path(dfo, f'orthologs/{gene}.fas'))
    paralogs_path = str(Path(dfo, f'paralogs/{gene}_paralogs.fas'))

    # parse orthologs for a given gene into a dictionary
    orthologs = SeqIO.to_dict(SeqIO.parse(orthologs_path, "fasta"))

    # parse paralogs into a dictionary (if they exist); else make an empty dict
    if os.path.isfile(paralogs_path):
        paralogs = SeqIO.to_dict(SeqIO.parse(paralogs_path, "fasta"))
    else:
        paralogs = {}

    # collect all candidate sequences for a given gen from the fisher.py result
    seq_dict = collect_seqs(gene)

    with open(table, 'r') as infile:
        for line in infile:
            branch_label, _, status = line.strip().split('\t')
            _, unique_id = branch_label.split('@')

            # Sequences already present in database
            if branch_label.count('_') == 1:
                record = seq_dict[unique_id]
                
                # Originally an ortholog
                if '..' not in branch_label:
                    # delete sequence all together
                    if status == 'd':
                        del orthologs[unique_id]
                    # change status from paralog to ortholog
                    elif status == 'p':
                        # prepare paralog name
                        pname = paralog_name(unique_id, paralogs.keys())
                        paralogs[pname] = record
                        # delete sequence from orthologs
                        del orthologs[unique_id]
                    else:
                        raise UnknownStatusError(unique_id, table)

                # Originally a paralog
                else:
                    unique_id, _ = unique_id.split('..')
                    # delete sequence all together
                    if status == 'd':
                        del paralogs[unique_id]
                    elif status == 'o':
                        record = paralogs[unique_id]
                        # for cases when ortholog has not survived trimming
                        if unique_id in orthologs:
                            # prepare paralog name
                            pname = paralog_name(unique_id, paralogs.keys())
                            # chance status from paralog to ortholog
                            paralogs[pname] = orthologs[unique_id]
                            # record in orthologs will be replaced
                            # with a new one in the next step
                        # add sequence to orthologs 
                        orthologs[unique_id] = record
                        # delete sequence from paralogs
                        del paralogs[unique_id]
                    else:
                        raise UnknownStatusError(unique_id, table)

            # New sequences
            elif branch_label.count('_') == 3:
                _, _, quality, _ = branch_label.split('_')
                qname = f'{unique_id}_{quality}'
                record = seq_dict[qname]
                if status == 'o':
                    # add to orthologs
                    orthologs[unique_id] = record
                elif status == 'p':
                    # add to paralogs
                    pname = paralog_name(unique_id, paralogs.keys())
                    paralogs[pname] = record

            
    return orthologs, paralogs


def add_to_meta(abbrev):
    """"Transfer input metadata for a given organism
    to dataset metadata.
    input:  short name of organism
    return: None
    """
    new_row = {'Unique ID': abbrev,
               'Long Name': input_info[abbrev]['full_name'],
               'Higher Taxonomy': input_info[abbrev]['tax'].replace('*', ''),
               'Lower Taxonomy': input_info[abbrev]['subtax'].replace('*', ''),
               'Data Type': input_info[abbrev]['data_type'],
               'Source': input_info[abbrev]['notes']
               }

    df = pd.read_csv(metadata, delimiter='\t')
    df.index = df['Unique ID']

    df = df.append(new_row, ignore_index=True)
    df = df.sort_values(by=['Higher Taxonomy', 'Lower Taxonomy', 'Unique ID'])
    df.to_csv(str(Path(dfo, 'metadata.tsv')), sep='\t', index=False)


def new_database(table):
    """Make changes in the PhyloFisher dataset according to information from
     parsed single gene tree."""

    gene = table.split('/')[-1].split('_parsed')[0]
    # Orthologs for a given gene
    orthologs_path = str(Path(dfo, f'orthologs/{gene}.fas'))
    # Paralogs for a given gene
    paralogs_path = str(Path(dfo, f'paralogs/{gene}_paralogs.fas'))

    orthologs, paralogs = parse_table(table)

    with open(orthologs_path, 'w') as res:
        # make changes in orthologs
        for name, record in orthologs.items():
            if name not in meta_orgs and name not in to_exclude:
                meta_orgs.add(name)
                add_to_meta(name)
            if name.split('..')[0] not in to_exclude:
                res.write(f'>{name}\n{record.seq}\n')

    with open(paralogs_path, 'w') as res:
        # make changes in paralogs
        for name, record in paralogs.items():
            if name.split('..')[0] not in meta_orgs and name.split('..')[0] not in to_exclude:
                meta_orgs.add(name.split('.')[0])
                add_to_meta(name.split('.')[0])
            if name.split('..')[0] not in to_exclude:
                res.write(f'>{name}\n{record.seq}\n')


def rebuild_db():
    """
    
    :return: 
    """
    try:
        shutil.rmtree(f'{dfo}/profiles')
    except FileNotFoundError:
        pass

    try:
        shutil.rmtree(f'{dfo}/datasetdb')
    except FileNotFoundError:
        pass
    
    cwd = os.getcwd()
    os.chdir(dfo)
    args.rename = None
    build_database.main(args, args.threads, True, 0.1)
    os.chdir(cwd)


def cp_proteomes():
    with open(input_metadata, 'r') as infile:
        infile.readline()
        file_dict = {}
        for line in infile:
            s_line = line.strip().split('\t')
            file_dict[s_line[2]] = os.path.abspath(os.path.join(s_line[0], s_line[1]))

    os.mkdir('tmp')
    os.chdir('tmp')
    for key in file_dict.keys():
        with tarfile.open(f'{key}.faa.tar.gz', "x:gz") as tar:
            tar.add(file_dict[key], arcname=f'{key}.faa')
        shutil.copy(f'{key}.faa.tar.gz', f'{dfo}/proteomes/{key}.faa.tar.gz')

    os.chdir('..')
    shutil.rmtree('tmp')


def main():
    """
    Main function. Run new_database on all parsed trees (tsv files)
    :return:
    """
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

    tools.backup(dfo)
    
    if args.to_exclude:
        to_exclude = taxa_to_exclude()
    else:
        to_exclude = []
        
    input_metadata = os.path.abspath(config['PATHS']['input_file'])
    metadata = str(Path(dfo, 'metadata.tsv'))
    meta_orgs = dataset_orgs()
    input_info = parse_input(input_metadata, to_exclude)
    main()
