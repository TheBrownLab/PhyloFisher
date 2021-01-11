#!/usr/bin/env python
import glob
import math
import os
import random
import shutil
import tarfile
import textwrap

import numpy
import pandas as pd

from phylofisher import help_formatter


def make_tmp_dir(base):
    """
    Makes tmpo directory. If it exists the old tmp directory is first removed
    """
    tmp_dir = f'{base}/tmp'
    try:
        os.mkdir(tmp_dir)
    except FileExistsError:
        shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir)


def get_num_of_replicates(x):
    """
    Uses 0.95=1-(1-x/100)^n to determine the number of replicates where x is the percentage of genes samples and
    n is the number of replicates
    """
    bottom_eq = 1 - (x / 100)
    confidence = 1 - args.confidence_interval
    n = numpy.log(confidence) / numpy.log(bottom_eq)
    n = math.ceil(n)
    return n


def subsample(percent, genes, output):
    """
    Creates gene subsets and super matirices of the subsampled gene sets.
    """
    replicate_num = get_num_of_replicates(percent)
    # Creates gene subsample sets with the number of replicates calculated in get_num_of_replicates
    for i in range(replicate_num):
        size = len(genes) * percent * 0.01
        sample = random.sample(genes, int(size))
        make_tmp_dir(output)
        # copy gene files into a tmp directory
        for file in sample:
            shutil.copy(file, f'{output}/tmp')
        # Runs forge to create super matrix of subsampled gene sets
        os.system(f'matrix_constructor.py '
                  f'-i {output}/tmp '
                  f'-o {output}/replicate_{i + 1} '
                  f'-if {args.in_format.lower()} '
                  f'-of {args.out_format.lower()} '
                  f'-c')
        shutil.rmtree(f'{output}/tmp')


def make_csv():
    perc_dict = dict()
    for x in sampling_percents:
        perc_dict[str(x)] = dict()
    for root, dirs, files in os.walk(args.output):
        for file in files:
            if file == 'indices.tsv':
                path = os.path.join(root, file)
                _, perc, rep, _ = path.split('/')
                perc, _ = perc.split('_')
                _, rep = rep.split('_')

                with open(path, 'r') as infile:
                    perc_dict[str(perc)][rep] = []
                    infile.readline()
                    for line in infile:
                        line = line.strip()
                        gene, _, _ = line.split('\t')
                        perc_dict[str(perc)][rep].append(gene)

    for percent in perc_dict.keys():
        data = perc_dict[percent]
        data = {int(k): v for k, v in data.items()}
        df = pd.DataFrame(data).transpose()
        df = df.sort_index().transpose()
        df = df.add_prefix('Replicate ')

        df.to_csv(f'{args.output}/{percent}_Percent.tsv', index=False, sep='\t')


def clean_up():
    """
    Copies matricies in to root output directory and tar's "percent" directories
    """
    for root, dirs, files in os.walk(args.output):
        for file in files:
            if file.endswith('.fas') or file.endswith('.phy'):
                path = os.path.join(root, file)
                ext = path.split('.')[-1]
                _, perc, rep, _ = path.split('/')
                perc, _ = perc.split('_')
                _, rep = rep.split('_')

                new_name = f'{args.output}/{perc}rep{rep}.{ext}'
                shutil.copy(path, new_name)

    os.chdir(args.output)
    for root, dirs, files in os.walk('.'):
        for dir in dirs:
            if dir.endswith('Percent'):
                with tarfile.open(f'{dir}.tar.gz', mode='x:gz') as my_tar:
                    my_tar.add(dir)
                shutil.rmtree(dir)


if __name__ == '__main__':
    description = 'Constructs super-matrices from randomly sampled genes.'
    parser, optional, required = help_formatter.initialize_argparse(name='random_resampler.py',
                                                                    desc=description,
                                                                    usage='random_resampler.py '
                                                                          '[OPTIONS] -i /path/to/input/')

    # Optional Arguments
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                              Format of the input single gene alignments.
                              Options: fasta, phylip (names truncated at 10 characters), 
                              phylip-relaxed (names are not truncated), or nexus.
                              Default: phylip-relaxed
                              """))
    optional.add_argument('-of', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                              Desired format of the output steps.
                              Options: fasta, nexus, phylip (names truncated at 10 characters), 
                              or phylip-relaxed (names are not truncated)
                              Default: phylip-relaxed
                              """))
    optional.add_argument('-ci', '--confidence_interval', metavar='0.N', type=float, default=0.95,
                          help=textwrap.dedent("""\
                              Confidence interval to use to calculate the number of replicates required.
                              Default: 0.95%%
                              """))
    optional.add_argument('-ps', '--percent_sampling', metavar='N', type=int, default=20,
                          help=textwrap.dedent("""\
                              Percent sampling step size.
                              Default: 20%% (e.g. 5%%, 10%%, 20%%, 25%%)
                              The default 20%% sampling results in a sampling series of 20%%, 40%%, 60%%, and 80%%.
                              """))

    in_help = 'Path to input directory containing gene files in FASTA format.'
    args = help_formatter.get_args(parser, optional, required, in_help=in_help)

    args.output = f'{args.output}_{args.percent_sampling}ps_{args.confidence_interval}ci'

    os.mkdir(args.output)
    # creates a list of percents to sample gene set at
    sampling_percents = list(range(0, 100, args.percent_sampling))[1:]
    # gets list of all gene files in the input directory
    files = [file for file in glob.glob(f'{args.input}/{args.prefix}*{args.suffix}')]

    # iterated through list of sampling percentages and creates a directory to store the outputs
    for x in sampling_percents:
        new_dir = f'{args.output}/{str(x)}_Percent'
        os.mkdir(new_dir)
        subsample(x, files, new_dir)

    make_csv()
    clean_up()
