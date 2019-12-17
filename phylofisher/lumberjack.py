#!/usr/bin/env python
import os
import glob
import string
import random
import argparse
from Bio import SeqIO
import configparser
from pathlib import Path
from collections import defaultdict


def parse_metadata():
    orgs = set()
    with open(metadata) as md:
        next(md)
        for line_ in md:
            sline = line_.split('\t')
            abbrev = sline[0].strip()
            orgs.add(abbrev)
    return orgs


def parse_input(multi_input):
    input_info = defaultdict(dict)
    for line in open(multi_input):
        sline = line.split('\t')
        abbrev = sline[2].strip()
        group = sline[3].strip()
        full_name = sline[6].strip()
        subtax = sline[4]
        col = sline[7]
        data_type = sline[8]
        notes = sline[9].strip()
        input_info[abbrev]['tax'] = group
        input_info[abbrev]['full_name'] = full_name
        input_info[abbrev]['subtax'] = subtax
        input_info[abbrev]['col'] = col
        input_info[abbrev]['data_type'] = data_type
        input_info[abbrev]['notes'] = notes
    return input_info



def parse_fasta(gene):
    #all seqs from fisher.py
    res = {}
    for record in SeqIO.parse(f'fasta/{gene}.fas', 'fasta'):
        if record.name.count('_') == 3:
            abbrev, _, _, quality = record.name.split('_')
            qname = f'{abbrev}_{quality}'
            res[qname] = record
        else:
            res[record.name] = record
    return res


def id_generator(size=5, chars=string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def paralog_name(abbrev, keys):
    id_ = id_generator()
    pname = f'{abbrev}..p{id_}'
    if pname not in keys:
        return pname
    else:
        paralog_name(abbrev, keys)


def parse_table(table):
    gene = table.split('/')[-1].split('_')[0]
    orthologs_path = str(Path(dfo, f'orthologs/{gene}.fas'))
    paralogs_path = str(Path(dfo, f'paralogs/{gene}_paralogs.fas'))
    orthologs = SeqIO.to_dict(SeqIO.parse(orthologs_path, "fasta"))
    if os.path.isfile(paralogs_path):
        paralogs = SeqIO.to_dict(SeqIO.parse(paralogs_path, "fasta"))
    else:
        paralogs = {}

    seq_dict = parse_fasta(gene)

    for line in open(table):
        # for orthologs from dataset
        tree_name, tax, status = line.split('\t')
        status = status.strip()
        abbrev = tree_name.split('@')[-1]
        if tree_name.count('_') != 3 and '..' not in abbrev:
            record = seq_dict[abbrev]
            if status == 'd':
                del orthologs[abbrev]
            elif status == 'p':
                pname = paralog_name(abbrev, paralogs.keys())
                paralogs[pname] = record
                del orthologs[abbrev]

    for line in open(table):
        #for new sequences
        tree_name, tax, status = line.split('\t')
        status = status.strip()
        abbrev = tree_name.split('@')[-1]
        if tree_name.count('_') == 3:
            quality = tree_name.split('_')[2]
            qname = f'{abbrev}_{quality}'
            record = seq_dict[qname]
            if status == 'o':
                orthologs[abbrev] = record
            elif status == 'p':
                pname = paralog_name(abbrev, paralogs.keys())
                paralogs[pname] = record

    for line in open(table):
        # for paralogs from dataset
        tree_name, tax, status = line.split('\t')
        status = status.strip()
        name = tree_name.split('@')[-1]
        abbrev = name.split('.')[0]
        if '..' in name:
            record = paralogs[name]
            if status == 'o':
                orthologs[abbrev] = record
                del paralogs[name]
            elif status == 'd':
                del paralogs[name]

    return orthologs, paralogs


def add_to_meta(abbrev):
    with open(metadata, 'a') as res:
        full = input_info[abbrev]['full_name']
        tax = input_info[abbrev]['tax']
        subtax = input_info[abbrev]['subtax']
        col = input_info[abbrev]['col']
        data_type = input_info[abbrev]['data_type']
        notes = input_info[abbrev]['notes']
        res.write(f'{abbrev}\t{full}\t{tax}\t{subtax}\t{col}\t{data_type}\t{notes}\n')


def new_database(table):
    gene = table.split('/')[-1].split('_')[0]
    print(str(Path(dfo, f'orthologs/{gene}.fas')))
    orthologs_path = str(Path(dfo, f'orthologs/{gene}.fas')) #Orthologs for a given gene
    paralogs_path = str(Path(dfo, f'paralogs/{gene}_paralogs.fas')) #Paralogs for a given gene

    orthologs, paralogs = parse_table(table)

    with open(orthologs_path, 'w') as res:
        for name, record in orthologs.items():
            if name not in meta_orgs:
                meta_orgs.add(name)
                add_to_meta(name)
            res.write(f'>{name}\n{record.seq}\n')

    with open(paralogs_path, 'w') as res:
        for name, record in paralogs.items():
            if name.split('.')[0] not in meta_orgs:
                meta_orgs.add(name.split('.')[0])
                add_to_meta(name.split('.')[0])
            res.write(f'>{name}\n{record.seq}\n')


def main():
    for table in glob.glob(f'{args.input_folder}/*.tsv'):
        new_database(table)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='lumberjack', usage="bla bla")
    parser.add_argument('-i', '--input_folder', required=True, help='folder with tsv files')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
    multi_input = os.path.abspath(config['PATHS']['input_file'])
    metadata = str(Path(dfo, 'metadata.tsv'))
    meta_orgs = parse_metadata()
    input_info = parse_input(multi_input)
    main()
