#!/usr/bin/env python
import os
import argparse
import glob
import textwrap
from collections import defaultdict
from Bio import SeqIO
import statistics
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from shutil import copyfile
import configparser

from phylofisher import fisher

"""
How are we going to get trimmed alignments?
"""


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


# Do we want to prepare trimmed alignments?
def prepare_alignments():
    pass


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


def extract_genes(s, dataset, gene_n, out_fold):
    genes = list(s.index[:gene_n])
    os.mkdir(out_fold)
    ortho_fold = os.path.join(dataset, 'orthologs')
    for gene in genes:
        g_name = f'{gene}.fas'
        copyfile(f'{os.path.join(ortho_fold, g_name)}', f'{os.path.join(out_fold, g_name)}')


def make_plot(s, plot_name, bin_size):
    means = []
    missing = []
    lengths = []
    matplotlib.use('pdf')
    for i, values in enumerate(s):
        if not means:
            missing.append(values[0])
            means.append(values[0])
            lengths.append(values[1])
        else:
            missing.append(values[0])
            means.append((sum(missing) / len(missing)))
            lengths.append((lengths[i - 1] + values[1]))

    labels_ = np.arange(0, len(lengths), bin_size)
    plt.plot(lengths, means)
    plt.xticks(lengths[::bin_size], labels=labels_)
    plt.ylabel('missing data [%]', labelpad=15)
    plt.xlabel('gene number', labelpad=15)
    plt.tight_layout()
    plt.savefig(f'{plot_name}.pdf')


def main(taxa, dataset, suffix, plot_name, bin_size, gene_n, out_fold):
    target_names = collect_names(taxa)
    metadata = os.path.join(dataset, 'metadata.tsv')
    tax_groups = parse_metadata(metadata)
    genes = {}
    for file in glob.glob(f'*{suffix}'):
        core = file.split('.')[0]
        genes[core] = missing_data(target_names, file, tax_groups)

    s = pd.Series(genes).sort_values()

    if plot_name:
        make_plot(s, plot_name, bin_size)

    if gene_n:
        extract_genes(s, dataset, gene_n, out_fold)


if __name__ == '__main__':
    formatter = lambda prog: fisher.myHelpFormatter(prog, max_help_position=100)

    parser = argparse.ArgumentParser(prog='missing_data.py',
                                     # TODO: Description
                                     description='some description',
                                     usage='missing_data.py [OPTIONS] -t <taxa> -d <dataset> -n <gene_number>',
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent("""\
                                     additional information:
                                        stuff
                                        """))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    # Required Arguments
    required.add_argument('-t', '--taxa', required=True, metavar='<taxa>',
                          # TODO: Description
                          help=textwrap.dedent("""\
                          some description  
                          """))
    required.add_argument('-d', '--dataset', required=True, metavar='<dataset>',
                          help=textwrap.dedent("""\
                          Path to input directory
                          """))
    required.add_argument('-n', '--gene_number', type=int, required=True, metavar='<N>',
                          help=textwrap.dedent("""\
                          Number of genes for analysis
                          """))

    # Optional Arguments
    optional.add_argument('-f', '--output_folder', type=str, default='output', metavar='<out_dir>',
                          help=textwrap.dedent("""\
                          Path to output directory
                          """))
    optional.add_argument('-s', '--suffix', type=str, metavar='<suff>',
                          help=textwrap.dedent("""\
                          Suffix of input files
                          Default: NONE
                          Example: path/to/input/*.suffix
                          """))
    optional.add_argument('-c', '--use_config', action='store_true',
                          help=textwrap.dedent("""\
                          Use config file for metadata
                          """))
    optional.add_argument('-p', '--plot',
                          help=textwrap.dedent("""\
                          Plot missing data statistics. 
                          By default statistics will NOT be plotted.
                          """))
    optional.add_argument('-b', '--bin_size', type=int, default=5, metavar='<N>',
                          help=textwrap.dedent("""\
                          Bin size for plot
                          Default: 5
                          """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                          Show this help message and exit.
                          """))

    parser._action_groups.append(optional)
    args = parser.parse_args()

    if args.use_config is True:
        config = configparser.ConfigParser()
        config.read('config.ini')
        args.dataset = str(Path(config['PATHS']['dataset_folder']).resolve())

    main(args.taxa, args.dataset, args.suffix,
         args.plot, args.bin_size, args.gene_number, args.output_folder)
