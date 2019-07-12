#!/usr/bin/env python
import os
from glob import glob
from Bio import SeqIO
from collections import defaultdict
import argparse
import csv


def parse_names(input_folder):
    name_set = set()
    if args.sufix:
        files = sorted(glob(f'{input_folder}/*{args.sufix}'))
    else:
        files = sorted(glob(f'{input_folder}/*'))
    for file in files:
        with open(file) as f:
            for line in f:
                if line.startswith('>'):
                    fname = line.split()[0][1:]
                    name = fname.split('_')[0]
                    name_set.add(name)
    return files, sorted(list(name_set))


def stats(total_len):
    with open('forge_stats.tsv', 'w') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['org', 'missing[%]'])
        missing = []
        for record in SeqIO.parse(args.output, 'fasta'):
            missing.append((record.name, (record.seq.count('-')/total_len)*100))
        for org_missing in sorted(missing, key=lambda x: x[1], reverse=True):
            tsv_writer.writerow(list(org_missing))

def main():
    input_folder = os.path.basename((args.input_folder.strip('/')))
    files, orgs = parse_names(input_folder)
    total_len = 0
    res_dict = defaultdict(str)
    for file in files:
        gene = os.path.basename(file).split('.')[0]
        length = 0
        seq_dict = {}
        for record in SeqIO.parse(file, 'fasta'):
            length = len(record.seq)
            seq_dict[record.id.split('_')[0]] = str(record.seq)
        total_len += length
        print(gene, total_len)
        for org in orgs:
            if org in seq_dict:
                res_dict[org] += seq_dict[org]
            else:
                res_dict[org] += ('-' * length)

    with open(args.output, "w") as res:
        for org, seq in res_dict.items():
            res.write(f'>{org}\n{seq}\n')
    stats(total_len)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='informant', usage="informant.py -i input_folder [OPTIONS]")
    parser.add_argument('-i', '--input_folder', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-s', '--sufix')
    args = parser.parse_args()
    main()