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


def parse_metadata(metadata_table):
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
        full_name = sline[5].strip()
        input_info[abbrev]['tax'] = group
        input_info[abbrev]['full_name'] = full_name
    return input_info


def check_mistakes():
    pass


def parse_fasta(gene):
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
    gene_meta_p = str(Path(dfo, f'genes/{gene}.fas'))
    para_meta_p = str(Path(dfo, f'paralogs/{gene}_paralogs'))
    gene_meta = SeqIO.to_dict(SeqIO.parse(gene_meta_p, "fasta"))
    if os.path.isfile(para_meta_p):
        para_meta = SeqIO.to_dict(SeqIO.parse(para_meta_p, "fasta"))
    else:
        para_meta = {}

    seq_dict = parse_fasta(gene)

    for line in open(table):
        tree_name, tax, status = line.split('\t')
        status = status.strip()
        abbrev = tree_name.split('@')[-1]
        if tree_name.count('_') != 3 and '..' not in tree_name:
            record = seq_dict[abbrev]
            if status == 'd':
                del gene_meta[abbrev]
            elif status == 'p':
                pname = paralog_name(abbrev, para_meta.keys())
                para_meta[pname] = record
                del gene_meta[abbrev]

    for line in open(table):
        tree_name, tax, status = line.split('\t')
        status = status.strip()
        abbrev = tree_name.split('@')[-1]
        if tree_name.count('_') == 3:
            quality = tree_name.split('_')[2]
            qname = f'{abbrev}_{quality}'
            record = seq_dict[qname]
            if status == 'o':
                gene_meta[abbrev] = record
            elif status == 'p':
                pname = paralog_name(abbrev, para_meta.keys())
                para_meta[pname] = record

    for line in open(table):
        tree_name, tax, status = line.split('\t')
        status = status.strip()
        name = tree_name.split('@')[-1]
        abbrev = name.split('.')[0]
        if '..' in tree_name:
            record = seq_dict[name]
            if status == 'o':
                gene_meta[abbrev] = record
            elif status == 'd':
                del gene_meta[name]

    return gene_meta, para_meta


def add_to_meta(abbrev):
    with open(metadata, 'a') as res:
        full = input_info[abbrev]['full_name']
        tax = input_info[abbrev]['tax']
        res.write(f'{abbrev}\t{full}\t{tax}\tx\tx\tx\n')


def new_database(table):
    gene = table.split('/')[-1].split('_')[0]
    print(str(Path(dfo, f'genes/{gene}.fas')))
    gene_meta_p = str(Path(dfo, f'genes/{gene}.fas'))
    para_meta_p = str(Path(dfo, f'paralogs/{gene}_paralogs'))

    gene_meta, para_meta = parse_table(table)

    with open(gene_meta_p, 'w') as res:
        for name, record in gene_meta.items():
            if name not in meta_orgs:
                meta_orgs.add(name)
                add_to_meta(name)
            res.write(f'>{name}\n{record.seq}\n')

    with open(para_meta_p, 'w') as res:
        for name, record in para_meta.items():
            if name.split('.')[0] not in meta_orgs:
                meta_orgs.add(name.split('.')[0])
                add_to_meta(name.split('.')[0])
            res.write(f'>{name}\n{record.seq}\n')


def main():
    for table in glob.glob(f'{args.input_folder}/*.tsv'):
        new_database(table)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='lumberjack', usage="bla bla")
    parser.add_argument('-i', '--input_folder', required=True)
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
    multi_input = os.path.abspath(config['PATHS']['input_file'])
    metadata = str(Path(dfo, 'metadata.tsv'))
    meta_orgs = parse_metadata(str(Path(dfo, 'metadata.tsv')))
    input_info = parse_input(multi_input)
    main()
