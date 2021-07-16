#!/usr/bin/env python
import configparser
import csv
import glob
import os
import shutil
import subprocess
import tarfile
import textwrap
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from phylofisher import help_formatter

SNAKEFILE_PATH = f'{os.path.dirname(os.path.realpath(__file__))}/sgt_constructor.smk'

def bash(cmd):
    """
    Function to run bash commands in a shell
    :param cmd:
    :return:
    """
    subprocess.run(cmd, shell=True, executable='/bin/bash')  # ,


def read_full_proteins(core):
    full_prots = {}
    for record in SeqIO.parse(f'{args.output}/prequal/{core}.aa.filtered', 'fasta'):
        full_prots[record.name] = record.seq
    return full_prots


def good_length(trimmed_aln, threshold):
    core = trimmed_aln.split('.')[0]
    full_proteins = read_full_proteins(core)
    original_name = f'{core}.length_filtered'
    length = None
    with open(original_name, 'w') as res:
        for record in SeqIO.parse(trimmed_aln, 'fasta'):
            if length is None:
                length = len(record.seq)
            coverage = len(str(record.seq).replace('-', '').replace('X', '')) / len(record.seq)
            if coverage > threshold:
                res.write(f'>{record.description}\n{full_proteins[record.name]}\n')
            else:
                print('deleted:', record.name, coverage)


def get_genes():
    return [file.split('.')[0] for file in os.listdir(args.input)]


def get_outfiles():
    if args.no_trees:
        outfiles = [f'{args.output}/trimal/{gene}.final' for gene in get_genes()]
    else:
        outfiles = [f'{args.output}-local.tar.gz']

    return ' '.join(outfiles)
    

def make_config():
    config_frags = [
        f'out_dir={args.output}',
        f'in_dir={args.input}',
        f'genes={",".join(get_genes())}',
        f'metadata={args.metadata}',
        f'tree_colors={args.color_conf}',
        f'input_metadata={args.input_metadata}',
        f'trees_only={args.trees_only}',
        f'no_trees={args.no_trees}'
    ]

    return ' '.join(config_frags)


def run_snakemake():
    smk_frags = [
        f'snakemake',
        f'-s {SNAKEFILE_PATH}',
        f'--config {make_config()}',
        f'--cores {args.threads}',
        f'--rerun-incomplete',
        f'--keep-going',
        f'--nolock'
    ]
    smk_cmd = ' '.join(smk_frags)
    smk_cmd += ' ' + get_outfiles()
    bash(smk_cmd)


def compress_output():
    """

    :return:
    """
    os.chdir(cwd)

    source_dir = f'tmp/{args.output}-local'

    shutil.copytree(f'{args.output}/trees', f'{source_dir}/trees')
    shutil.copy(args.metadata, f'{source_dir}/{os.path.basename(args.metadata)}')
    shutil.copy(args.input_metadata, f'{source_dir}/{os.path.basename(args.input_metadata)}')
    shutil.copy(args.color_conf, f'{source_dir}/{os.path.basename(args.color_conf)}')

    with tarfile.open(f'{args.output}-local.tar.gz', "x:gz") as tar:
        tar.add(source_dir, arcname=f'{args.output.split("/")[-1]}-local')

    shutil.rmtree('tmp')


if __name__ == '__main__':
    description = 'Aligns, trims, and builds single gene trees from unaligned gene files.'
    usage = 'sgt_constructor.py -i path/to/input/ [OPTIONS]'
    parser, optional, required = help_formatter.initialize_argparse(name='sgt_constructor.py',
                                                                    desc=description,
                                                                    usage=usage)

    # Optional Arguments
    optional.add_argument('-t', '--threads', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Desired number of threads to be utilized.
                          Default: 1
                          """))
    optional.add_argument('--no_trees', action='store_true',
                          help=textwrap.dedent("""\
                          Do NOT build single gene trees.
                          Length filtration and trimmming only.
                          """))
    optional.add_argument('--trees_only', action='store_true',
                          help=textwrap.dedent("""\
                          Only build single gene trees.
                          No length filtration and trimming.
                          """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Format of the input files.
                          Options: fasta, phylip (names truncated at 10 characters), 
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=True)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    args.metadata = str(os.path.join(dfo, 'metadata.tsv'))
    args.color_conf = str(os.path.abspath(config['PATHS']['color_conf']))
    args.input_metadata = str(os.path.abspath(config['PATHS']['input_file']))
    args.input = str(os.path.abspath(args.input[:-1])) if args.input.endswith('/') else str(os.path.abspath(args.input))
    args.output = str(os.path.abspath(args.output[:-1])) if args.output.endswith('/') else str(os.path.abspath(args.output))

    out_dict = {'fasta': 'fas',
                'phylip': 'phy',
                'phylip-relaxed': 'phy',
                'nexus': 'nex'}

    run_snakemake()