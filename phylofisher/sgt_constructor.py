#!/usr/bin/env python
import configparser
import os
import shutil
import subprocess
import textwrap
from pathlib import Path
from collections import OrderedDict
from phylofisher import help_formatter

SNAKEFILE_PATH = f'{os.path.dirname(os.path.realpath(__file__))}/sgt_constructor.smk'


def get_genes(length_filter):
    '''
    Get list of genes from input dir

    :return: list of genes
    :rtype: list
    '''
    if length_filter or args.trees_only:
        ret = [file.split('.')[0] for file in os.listdir(args.input)]
    
    else:
        ret = []
        len_filt_bmge_dir = f'{args.output}/length_filtration/bmge'
        bmge_out_files = [file for file in os.listdir(len_filt_bmge_dir) if file.endswith('.bmge')]
        for bmge_out_file in bmge_out_files:
            with open(f'{len_filt_bmge_dir}/{bmge_out_file}', 'r') as infile:
                line = infile.readline()
                if line == '':
                    pass
                else:
                    ret.append(bmge_out_file.split('.bmge')[0])

    return ret
    

def make_config(length_filter):
    '''
    Make config list to be passed to 

    :return: snakemake config
    :rtype: list
    '''
    ret = [
        f'out_dir={args.output}',
        f'in_dir={args.input}',
        f'genes={",".join(get_genes(length_filter))}',
        f'trees_only={args.trees_only}',
        f'no_trees={args.no_trees}',
        f'database={args.database}',
        f'input_metadata={args.input_metadata}'
    ]

    return ' '.join(ret)


def get_output_files(length_filter):

    ret = []

    if length_filter:
        for gene in get_genes(length_filter):
            ret.append(f'{args.output}/length_filtration/bmge/{gene}.bmge')
    
    else:
        if args.no_trees:
            for gene in get_genes(length_filter):
                ret.append(f'{args.output}/trimal/{gene}.final')
        else:
            ret.append(f'{args.output}-local.tar.gz')
            
    return ' '.join(ret)


def run_snakemake(length_filter=False):
    '''
    Combine snakemake cmd frags and runs snakemake
    '''
    smk_frags = [
        f'snakemake',
        f'-s {SNAKEFILE_PATH}',
        f'--config {make_config(length_filter)}',
        f'--cores {args.threads}',
        f'--rerun-incomplete',
        f'--keep-going',
        f'--nolock',
        f'--use-conda',
    ]

    smk_frags.append(get_output_files(length_filter))
        
    smk_cmd = ' '.join(smk_frags)
    print(smk_cmd)
    subprocess.run(smk_cmd, shell=True, executable='/bin/bash')


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
    args.database = str(os.path.join(dfo, 'phylofisher.db'))
    args.input_metadata = str(os.path.abspath(config['PATHS']['input_file']))
    args.input = str(os.path.abspath(args.input[:-1])) if args.input.endswith('/') else str(os.path.abspath(args.input))
    args.output = str(os.path.abspath(args.output[:-1])) if args.output.endswith('/') else str(os.path.abspath(args.output))

    out_dict = {'fasta': 'fas',
                'phylip': 'phy',
                'phylip-relaxed': 'phy',
                'nexus': 'nex'}

    if not args.trees_only:
        run_snakemake(length_filter=True)

    run_snakemake(length_filter=False)