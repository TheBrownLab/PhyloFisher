#!/usr/bin/env python
import os, re
import shutil
import subprocess
import textwrap

import numpy as np
import pandas as pd

from phylofisher import help_formatter


def run_raxml_rtc():
    """
    :return:
    """
    if os.path.isdir(f'{args.output}/RAxML'):
        shutil.rmtree(f'{args.output}/RAxML')
    os.mkdir(f'{args.output}/RAxML')
    os.chdir(f'{args.output}/RAxML')
    genes = [file.split('.')[1] for file in os.listdir(args.input) if 'RAxML_bipartitions.' in file]
    genes = set(genes)

    for gene in genes:
        cmd = (f'raxmlHPC-PTHREADS-AVX2 '
               f'-t {args.input}/RAxML_bipartitions.{gene}{args.suffix} '
               f'-z {args.input}/RAxML_bootstrap.{gene}{args.suffix} '
               f'-n IC.{gene}{args.suffix} '
               f'-f i '
               f'-m PROTCATLGF '
               f'-T 2 '
               f'-s NA')

        subprocess.run(cmd, shell=True, executable='/bin/bash')

    os.chdir(cwd)


def get_rtc():
    """

    :return:
    """
    genes = [file.split('.')[1] for file in os.listdir(args.input) if 'RAxML_bipartitions.' in file]
    rtc_dict = {}
    for gene in genes:
        with open(f'{args.output}/RAxML/RAxML_info.IC.{gene}{args.suffix}', 'r') as infile:
            for line in infile:
                line = line.strip()
                if 'Tree certainty for this tree:' in line:
                    rtc = float(line.split(': ')[1])
                    rtc_dict[gene] = rtc

    return rtc_dict


def get_genes():
    """

    :return:
    """
    s = pd.Series(rtc_values)
    s = s.sort_values(ascending=False)
    # Get top 25%, 50%, and 75%
    s25 = s[s > np.percentile(s, 75)]
    s50 = s[s > np.percentile(s, 50)]
    s75 = s[s > np.percentile(s, 25)]
    # Get genes in each of the following pandas series
    genes_25 = s25.index
    genes_50 = s50.index
    genes_75 = s75.index

    return {'25': genes_25, '50': genes_50, '75': genes_75}


def build_matrices():
    gene_subsets = get_genes()
    gene_files = [os.path.join(args.input, file) for file in os.listdir(args.input)
                  if file.endswith(out_dict[args.in_format])]

    for key, subset in gene_subsets.items():
        if os.path.isdir(f'{args.output}/rtc{key}'):
            shutil.rmtree(f'{args.output}/rtc{key}')
        os.mkdir(f'{args.output}/rtc{key}')
    
        for gene_file in gene_files:
            for gene in subset:
                if gene in gene_file:
                    dst = f'{args.output}/rtc{key}/{os.path.basename(gene_file)}'
                    shutil.copy(gene_file, dst)

        cmd = (f'matrix_constructor.py '
               f'-i {args.output}/rtc{key} '
               f'-o {args.output}/matrix_constructor{key} '
               f'-if {args.in_format.lower()} '
               f'-of {args.out_format.lower()} '
               f'-c')
        subprocess.run(cmd, shell=True, executable='/bin/bash')


if __name__ == '__main__':
    description = ('Bins orthologs based on relative tree certainty score (RTC) calculated in RAxML \n'
                   'and constructs supermatrices from the bins.')
    parser, optional, required = help_formatter.initialize_argparse(name='rtc_binner.py',
                                                                    desc=description,
                                                                    usage='rtc_binner.py '
                                                                          '[OPTIONS] -i /path/to/input/')

    # Optional Arguments
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Format of the input alignments.
                          Options: fasta, nexus, phylip (names truncated at 10 characters),
                          or phylip-relaxed (names are not truncated)
                          Default: fasta
                          """))

    optional.add_argument('-of', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Desired format of the RTC binned alignments.
                          Options: fasta, nexus, phylip (names truncated at 10 characters), 
                          or phylip-relaxed (names are not truncated)
                          Default: fasta
                          """))

    optional.add_argument('-s', '--suffix', metavar='<suffix>', type=str, default='NONE',
                          help=textwrap.dedent("""\
                          Suffix of RAxML output files
                          Example: RAxML_bipartitions.GENE.<SUFFIX>, RAxML_bootstrap.GENE.<SUFFIX>
                          Default: NONE
                          """))

    in_help = ('Path to directory containing single gene trees built from only orthologs, \n'
               'corresponding bootstrap value files, and corresponding alignments.')
    args = help_formatter.get_args(parser, optional, required, in_help=in_help, pre_suf=False)

    args.output = os.path.abspath(args.output)
    args.input = os.path.abspath(args.input)

    out_dict = {'fasta': 'fas',
                'phylip': 'phy',
                'phylip-relaxed': 'phy',
                'nexus': 'nex'}

    if args.suffix and args.suffix[0] != '.':
        args.suffix = f'.{args.suffix}'

    if os.path.isdir(args.output):
        shutil.rmtree(args.output)
    os.mkdir(args.output)

    cwd = os.getcwd()

    run_raxml_rtc()
    rtc_values = get_rtc()
    build_matrices()
