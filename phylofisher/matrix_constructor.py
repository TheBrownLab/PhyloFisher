#!/usr/bin/env python
import os
import subprocess
import textwrap
from collections import defaultdict
from glob import glob
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofisher import help_formatter

SNAKEFILE_PATH = f'{os.path.dirname(os.path.realpath(__file__))}/matrix_constructor.smk'

def get_genes():
    '''
    Get input gene files

    :return: input gene files
    :rtype: list
    '''
    genes = [file.split(".")[0] for file in os.listdir(args.input)]

    return genes

def make_config():
    '''
    Make config list to be passed to 

    :return: snakemake config
    :rtype: list
    '''
    ret = [
        f'out_dir={args.output}',
        f'in_dir={args.input}',
        f'in_format={args.in_format}',
        f'out_format={args.out_format}',
        f'concatenation_only={args.concatenation_only}',
        f'genes={",".join(get_genes())}',
    ]

    return ' '.join(ret)


def get_output_files():

    ret = [
        f'{args.output}/matrix.fas',
        f'{args.output}/indices.tsv',
        f'{args.output}/matrix_constructor_stats.tsv'
    ]

    return ' '.join(ret)

def clean_up():
    '''
    Clean up intermediate files
    '''
    files_to_remove = glob(f'{args.output}/prequal/*.PP')
    files_to_remove += glob(f'{args.output}/divvier/*.fas')
    for file in files_to_remove:
        print(f'Removing {file}')
        os.remove(file)
    

def run_snakemake(length_filter=False):
    '''
    Combine snakemake cmd frags and runs snakemake
    '''
    smk_frags = [
        f'snakemake',
        f'-s {SNAKEFILE_PATH}',
        f'--config {make_config()}',
        f'--cores {args.threads}',
        f'--rerun-incomplete',
        f'--keep-going',
        f'--use-conda',
    ]

    smk_frags.append(get_output_files())
        
    smk_cmd = ' '.join(smk_frags)
    print(smk_cmd)
    subprocess.run(smk_cmd, shell=True, executable='/bin/bash', check=True)


if __name__ == '__main__':
    description = 'To trim align and concatenate orthologs into a super-matrix run:'
    parser, optional, required = help_formatter.initialize_argparse(name='matrix_constructor.py',
                                                                    desc=description,
                                                                    usage='matrix_constructor.py [OPTIONS] -i path/to/input/')

    # Optional Arguments
    optional.add_argument('-of', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Desired format of the output matrix.
                          Options: fasta, phylip (names truncated at 10 characters), 
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Format of the input files.
                          Options: fasta, phylip (names truncated at 10 characters), 
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))
    optional.add_argument('-c', '--concatenation_only', action='store_true', default=False,
                          help=textwrap.dedent("""\
                          Only concatenate alignments. Trimming and alignment are not performed automatically.
                          """))
    optional.add_argument('-t', '--threads', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Desired number of threads to be utilized.
                          Default: 1
                          """))
    optional.add_argument('--clean_up', action='store_true', default=False,
                          help=textwrap.dedent("""\
                          Clean up large intermediate files.
                          """))

    # Changes help descriptions from the default input and output help descriptions
    in_help = 'Path to prep_final_dataset_<M.D.Y>'
    args = help_formatter.get_args(parser, optional, required, in_help=in_help)

    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)

    run_snakemake()

    if args.clean_up:
        clean_up()

    