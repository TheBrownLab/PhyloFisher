#!/usr/bin/env python
import glob
import argparse
from Bio import SeqIO
import configparser
from pathlib import Path


def check_mistakes():
    pass


def parse_fasta(gene):
    res = {}
    for record in SeqIO.parse(f'fasta/{gene}.fas', 'fasta'):
        if record.name.count('_') == 3:
            abbrev, _, _, quality = record.name.split()
            qname = f'{abbrev}_{quality}'
            res[qname] = record
        else:
            res[record.name] = record
    return res


def parse_table(table):
    gene = table.split('_')[0]
    seq_dict = parse_fasta(gene)
    for line in open(table):
        tree_name, tax, status = line.split('\t')
        status = status.strip()
        abbrev = tree_name.split('@')[-1]
        if tree_name.count('_') == 3:
            quality = tree_name.split('_')[2]
            qname = f'{abbrev}_{quality}'
            record = seq_dict[qname]
        else:
            record = seq_dict[abbrev]




# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Script for ortholog fishing.', usage="fisher [OPTIONS]")
#     parser.add_argument('-i', '--input_folder', required=True)
#     parser.add_argument('-o' '--output_folder', required=True)
#     parser.add_argument('--add_to_database', action='store_true')
#     args = parser.parse_args()
#
#     config = configparser.ConfigParser()
#     config.read('config.ini')
#     dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
#     multi_input = os.path.abspath(config['PATHS']['input_file'])
#
#     meta_orgs = metada_orgs(str(Path(dfo, 'metadata.tsv')))