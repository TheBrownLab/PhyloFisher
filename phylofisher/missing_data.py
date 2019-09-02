#!/usr/bin/env python
import argparse
import glob
from collections import defaultdict
from Bio import SeqIO
import statistics
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path



def parse_metadata(metadata):
    groups = defaultdict(set)
    for line in open(metadata):
        if line.strip():
            if "Full Name" not in line:
                org, _, group, subtax, _, _, _ = line.split('\t')
                groups[group].add(org)
                groups[subtax].add(org)
    return groups

def collect_names(infile):
    names = []
    for line in open(infile):
        names.append(line.strip())
    return names


def parse_alignments(aln_file):
    res_dict = {}
    length = None
    for record in SeqIO.parse(aln_file, 'fasta'):
        res_dict[record.name] = str(record.seq)
        if not length:
            length = len(str(record.seq))
    return res_dict, length



def missing_data(target_names, alignments_file, tax_groups):
    ali_dict, length = parse_alignments(alignments_file)
    missing = []
    for name in target_names:
        if name in tax_groups:
            group = tax_groups[name]
        else:
            group = [name]

        for member in group:
            if member not in ali_dict:
                missing.append(1)
            else:
                seq = ali_dict[member]
                missing.append((seq.count('-') + seq.upper().count('X')) / len(seq))
    return (statistics.mean(missing), length)



def main(taxa, metadata, suffix, output):
    target_names = collect_names(taxa)
    tax_groups = parse_metadata(metadata)
    genes = {}
    for file in glob.glob(f'*{suffix}'):
        core = file.split('.')[0]
        genes[core] = missing_data(target_names, file, tax_groups)

    s = pd.Series(genes).sort_values()

    means = []
    missing = []
    lengths = []
    for i, values in enumerate(s):
        print(i, values)
        if not means:
            missing.append(values[0])
            means.append(values[0])
            lengths.append(values[1])
        else:
            missing.append(values[0])
            means.append((sum(missing)/len(missing)))
            lengths.append((lengths[i-1] + values[1]))

    labels_ = np.arange(0,len(lengths),5)
    plt.plot(lengths, means)
    plt.xticks(lengths[::5], labels=labels_)
    plt.ylabel('missing data [%]', labelpad=15)
    plt.xlabel('gene number', labelpad=15)
    plt.tight_layout()
    plt.savefig(f'{output}.pdf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for ortholog fishing.', usage="fisher.py [OPTIONS]")
    parser.add_argument('-t', '--taxa', required=True)
    parser.add_argument('-m', '--metadata')
    parser.add_argument('-s', '--suffix', type=str)
    parser.add_argument('-o', '--output')
    parser.add_argument('-c', '--use_config', action='store_true', help='use config file for metadata')
    args = parser.parse_args()

    if args.use_config is True:
        config = configparser.ConfigParser()
        config.read('config.ini')
        dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
        args.metadata = os.path.join(dfo, 'metadata.tsv')


    main(args.taxa, args.metadata, args.suffix, args.output)
