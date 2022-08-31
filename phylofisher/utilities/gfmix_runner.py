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
        f'iqtree={args.iqtree}',
        f'rootfile={args.rootfile}',
        f'basename={args.basename}',
    ]

    return ' '.join(ret)


def get_output_files():
    '''
    Returns expected output files

    :return: output files
    :rtype: list
    '''
    ret = []
    ret.append(f'{args.output}/{args.basename}.loglikelihood')
            
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
    description = 'Runs gfmix'
    parser, optional, required = help_formatter.initialize_argparse(name='gfmix_runner.py',
                                                                    desc=description,
                                                                    usage='gfmix_runner.py '
                                                                          '[OPTIONS] -t <tree_file> '
                                                                          '-m <matrix> -iq <iqtree> '
                                                                          '-r <rootfile> -b <basename>')

    # Required Arguments
    required.add_argument('-m', '--matrix', required=True, type=str, metavar='<matrix>',
                          help=textwrap.dedent("""\
                              Path to supermatrix file.
                              """))
    required.add_argument('-tr', '--tree', required=True, type=str, metavar='<tree>',
                          help=textwrap.dedent("""\
                              Path to tree.
                              """))
    required.add_argument('-iq', '--iqtree', required=True, type=str, metavar='<tree>',
                          help=textwrap.dedent("""\
                              Path to iqtree.
                              """))
    required.add_argument('-r', '--rootfile', required=True, type=str, metavar='<tree>',
                          help=textwrap.dedent("""\
                              Path root file. A txt file containing a list of taxa on one side of the root split.
                              Example:
                                taxonA
                                taxonB
                              """))
    required.add_argument('-b', '--basename', required=True, type=str, metavar='<tree>',
                          help=textwrap.dedent("""\
                              Output file basename.
                              """))

    # Optional Arguments
    optional.add_argument('-f', '--frequencies', metavar='<N>', type=int, default=60,
                          help=textwrap.dedent("""\
                              A file with the frequencies for each frequency class as rows. This file be obtained using MAMMaL.
                              If no file is given. MAMMaL will be run with the default options. (i.e. 60 frequency classes in the mixture model.)
                              Default: None
                            """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                              Input format of matrix
                              Options: fasta, nexus, phylip (names truncated at 10 characters), 
                              or phylip-relaxed (names are not truncated)
                              Default: fasta
                              """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)
    run_snakemake()
   
