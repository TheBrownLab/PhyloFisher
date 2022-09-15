#!/usr/bin/env python3

import os
import subprocess
import textwrap

from phylofisher import help_formatter

SNAKEFILE_PATH = f'{os.path.dirname(os.path.realpath(__file__))}/gfmix_mammal.smk'


def bash(cmd):
    '''
    Runs command ins bash shell

    :param cmd: command
    :type cmd: str
    '''
    p = subprocess.run(cmd, executable='/bin/bash', shell=True)



def make_config():
    '''
    Make config list to be passed to 

    :return: snakemake config
    :rtype: list
    '''
    ret = [
        f'out_dir={args.output}',
        f'in_format={args.in_format}',
        f'matrix={args.matrix}',
        f'tree={args.tree}',
        f'iqtree=None',
        f'rootfile=None',
        f'basename={args.basename}',
        f'rate_classes={args.rate_classes}',
    ]

    return ' '.join(ret)


def get_output_files():
    '''
    Returns expected output files

    :return: output files
    :rtype: list
    '''
    ret = []
    
    ret.append(f'{args.output}/{args.basename}.estimated-frequencies')
    ret.append(f'{args.output}/{args.basename}.esmodel.nex')
            
    return ' '.join(ret)


def run_snakemake():
    '''
    Combine snakemake cmd frags and runs snakemake
    '''

    
    smk_frags = [
        f'snakemake',
        f'-s {SNAKEFILE_PATH}',
        f'--config {make_config()}',
        f'--rerun-incomplete',
        f'--cores 1',
        f'--keep-going',
        f'--use-conda',
        f'--conda-prefix /tmp',
        f'--nolock'
    ]

    smk_frags.append(get_output_files())
        
    smk_cmd = ' '.join(smk_frags)
    print(smk_cmd)
    bash(smk_cmd)



if __name__ == "__main__":
    description = 'Prepares supermatrix and tree for MAMMaL analysis'
    parser, optional, required = help_formatter.initialize_argparse(name='mammal_modeler.py',
                                                                    desc=description,
                                                                    usage='mammal_modeler.py '
                                                                          '[OPTIONS] -t <tree_file> -s <matrix>')

    # Required Arguments
    required.add_argument('-m', '--matrix', required=True, type=str, metavar='<matrix>',
                          help=textwrap.dedent("""\
                              Path to supermatrix file.
                              """))
    required.add_argument('-tr', '--tree', type=str, metavar='<tree>',
                          help=textwrap.dedent("""\
                              Path to tree.
                              """))
    # Optional Arguments
    optional.add_argument('-c', '--rate_classes', metavar='<N>', type=int, default=60,
                          help=textwrap.dedent("""\
                              The number of frequency classes in the mixture model.
                              Options: 10, 20, 30, 40, 50, or 60
                              Default: 60
                            """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                              Input format of matrix
                              Options: fasta, nexus, phylip (names truncated at 10 characters), 
                              or phylip-relaxed (names are not truncated)
                              Default: fasta
                              """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)
    args.basename = os.path.basename(args.matrix).split('.')[0]

    run_snakemake()
