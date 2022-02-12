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

def bash(cmd):
    '''
    Function to run bash commands in a shell

    :param cmd: Command to run
    :type cmd: str
    '''
    subprocess.run(cmd, shell=True, executable='/bin/bash')  # ,


def get_genes():
    '''
    Get list of genes from input dir

    :return: list of genes
    :rtype: list
    '''
    return [file.split('.')[0] for file in os.listdir(args.input)]
    

def make_config():
    '''
    Make config list to be passed to 

    :return: snakemake config
    :rtype: list
    '''
    
    config_frags = [
        f'out_dir={args.output}',
        f'in_dir={args.input}',
        f'genes={",".join(get_genes())}',
        f'trees_only={args.trees_only}',
        f'no_trees={args.no_trees}',
    ]
    return ' '.join(config_frags)


def run_snakemake():
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
        f'--nolock'
    ]

    if args.dry_run:
        smk_frags.append('-n')
        
    smk_cmd = ' '.join(smk_frags)
    bash(smk_cmd)


def gather_and_compress_output():
    '''
    Copies relevant files to local output dir and compresses it
    '''
    os.mkdir(f'{args.output}/trees')
    
    local_out_dir = f'{args.output}-local'
    os.mkdir(local_out_dir)
    os.mkdir(f'{local_out_dir}/trees')
    
    to_copy = OrderedDict()
    to_copy[f'{args.color_conf}'] = f'{args.output}-local/tree_colors.tsv'
    to_copy[f'{args.metadata}'] = f'{args.output}-local/metadata.tsv'
    to_copy[f'{args.input_metadata}'] = f'{args.output}-local/input_metadata.tsv'

    genes = get_genes()
    for gene in genes:
        # if trees_only
        if args.trees_only:
            to_copy[f'{args.output}/{gene}.fas'] = f'{local_out_dir}/trees/{gene}.final'
            to_copy[f'{args.output}/raxml/RAxML_bipartitions.{gene}.tre'] = f'{local_out_dir}/trees/RAxML_bipartitions.{gene}.tre'
            to_copy[f'{args.output}/{gene}.fas'] = f'{local_out_dir}-local/trees/{gene}.final'
            to_copy[f'{args.output}/trees/RAxML_bipartitions.{gene}.tre'] = f'{local_out_dir}-local/trees/RAxML_bipartitions.{gene}.tre'
        
        # if pre processing and trees
        else:
            to_copy[f'{args.output}/length_filtration/bmge/{gene}.length_filtered'] = f'{local_out_dir}/trees/{gene}.trimmed'
            to_copy[f'{args.output}/trimal/{gene}.final'] = f'{local_out_dir}/trees/{gene}.final'
            to_copy[f'{args.output}/raxml/RAxML_bipartitions.{gene}.tre'] = f'{args.output}/trees/RAxML_bipartitions.{gene}.tre'
            to_copy[f'{args.output}/trees/{gene}.trimmed'] = f'{local_out_dir}/trees/{gene}.trimmed'
            to_copy[f'{args.output}/trees/{gene}.final'] = f'{local_out_dir}/trees/{gene}.final'          
            to_copy[f'{args.output}/trees/RAxML_bipartitions.{gene}.tre'] = f'{local_out_dir}/trees/RAxML_bipartitions.{gene}.tre'

    for k, v in to_copy.items():
        # nested ifs are to keep empty bmge files out for genes that didn't make it through the checkpoint
        if os.path.isfile(k):
            if os.stat(k).st_size > 0:
                shutil.copyfile(k, v) 

    local_out_dir_base = f'{args.output.split("/")[-1]}-local'
    tar_cmd = f'tar -czf {args.output}-local.tar.gz {local_out_dir_base}'
    subprocess.run(tar_cmd, shell=True, executable='/bin/bash')

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
    optional.add_argument('--dry_run', action='store_true',
                          help=textwrap.dedent("""\
                          Builds DAG and show jobs that would be performed.
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
    gather_and_compress_output()