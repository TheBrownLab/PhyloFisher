#!/usr/bin/env python
import os
import argparse
import configparser
from pathlib import Path
from Bio import SeqIO


def parse_genes(gene_file):
    to_exlude = set()
    with open(gene_file) as lines:
        next(lines)
        for line in lines:
            gene, _, _, sgt = line.split(',')
            sgt = sgt.strip()
            if sgt.lower() != 'yes':
                to_exlude.add(gene)
    return to_exlude


def parse_orgs(org_file):
    to_exlude = set()
    paralogs = set()
    with open(org_file) as lines:
        header = next(lines)
        if header.count(',') == 10:
            sgt_idx = 9
        else:
            sgt_idx = 6
        for line in lines:
            org = line.split(',')[0]
            sgt = line.split(',')[sgt_idx].strip()
            if sgt.lower() != 'yes':
                to_exlude.add(org)
            if header.count(',') == 10:
                if line.split(',')[sgt_idx + 1].strip().lower() == "yes":
                    paralogs.add(org)
    return to_exlude, paralogs


def fasta_filtr(file, o_to_ex, paralogs=None):
    with open(str(Path(args.output_directory, file)), 'w') as res:
        for record in SeqIO.parse(str(Path(args.input_directory, file)), 'fasta'):
            if record.name.split('_')[0] not in o_to_ex:
                res.write(f'>{record.name}\n{record.seq}\n')
        if paralogs:
            para_file = str(Path(dfo, f'paralogs/{file.split(".")[0]}_paralogs.fas'))
            if os.path.isfile(para_file):
                for record in SeqIO.parse(para_file, 'fasta'):
                    if record.name.split('.')[0] in paralogs:
                        res.write(f'>{record.name}\n{record.seq}\n')


def main():
    if args.orthologs:
        gene_file = 'orthologs_stats/genes_stats.csv'
        orgs_file = 'orthologs_stats/orgs_stats.csv'
    else:
        gene_file = str(Path(args.input_directory + '_stats', 'genes_stats.csv'))
        orgs_file = str(Path(args.input_directory + '_stats', 'orgs_stats.csv'))
    g_to_ex = parse_genes(gene_file)
    o_to_ex, paralogs = parse_orgs(orgs_file)
    filtered_genes = []
    for file in os.listdir(args.input_directory):
        if file.split('.')[0] not in g_to_ex:
            filtered_genes.append(file)

    for file in filtered_genes:
        if paralogs:
            fasta_filtr(file, o_to_ex, paralogs)
        else:
            fasta_filtr(file, o_to_ex)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for filtering orgs [and|or] genes',
                                     usage="fishing_net.py -i input_directory [OPTIONS]")
    parser.add_argument('-i', '--input_directory')
    parser.add_argument('-o', '--output_directory', required=True)
    parser.add_argument('--orthologs', action='store_true')
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
    if args.orthologs:
        args.input_directory = str(Path(dfo, 'orthologs'))
    os.mkdir(args.output_directory)
    main()
