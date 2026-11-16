#!/usr/bin/env python
import configparser
import os
import subprocess
import sys
import textwrap
from pathlib import Path
from phylofisher import help_formatter

SNAKEFILE_PATH = f'{os.path.dirname(os.path.realpath(__file__))}/sgt_constructor.smk'


def get_genes():
    '''
    Get list of genes from input dir

    :return: list of genes
    :rtype: list
    '''
    return [file.split('.')[0] for file in os.listdir(args.input)]


def make_config():
    '''
    Make config list to be passed to snakemake

    :return: snakemake config
    :rtype: str
    '''
    ret = [
        f'out_dir={args.output}',
        f'in_dir={args.input}',
        f'genes={",".join(get_genes())}',
        f'trees_only={args.trees_only}',
        f'no_trees={args.no_trees}',
        f'length_filter_only={args.length_filter_only}',
        f'tree_colors={args.color_conf}',
        f'metadata={args.metadata}',
        f'input_metadata={args.input_metadata}',
        f'threads_per_job={args.threads_per_job}'
    ]

    return ' '.join(ret)


def get_output_files():
    '''
    Get final output files for the workflow

    :return: space-separated list of output files
    :rtype: str
    '''
    if args.length_filter_only:
        return f'{args.output}-length_filtered.tar.gz'
    elif args.no_trees:
        # Aggregate target; genes emptied during length filtration are skipped
        return f'{args.output}/trimal/no_trees.done'
    else:
        # Full pipeline with trees
        return f'{args.output}-local.tar.gz'


def run_snakemake():
    '''
    Combine snakemake cmd frags and runs snakemake
    '''
    smk_frags = [
        f'snakemake',
        f'-s {SNAKEFILE_PATH}',
        f'--config {make_config()}',
        f'--rerun-incomplete',
        f'--keep-going',
        f'--nolock',
        f'--use-conda',
    ]

    if args.profile is not None:
        smk_frags.append(f'--profile {args.profile}')
    else:
        smk_frags.append(f'--cores {args.threads}')

    if args.dry_run:
        smk_frags.append('-n')

    smk_frags.append(get_output_files())

    smk_cmd = ' '.join(smk_frags)
    print(smk_cmd)
    try:
        subprocess.run(smk_cmd, shell=True, executable='/bin/bash', check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command '{e.cmd}' failed with return code {e.returncode}.", file=sys.stderr)
        sys.exit(e.returncode)


if __name__ == '__main__':
    description = 'Aligns, trims, and builds single gene trees from unaligned gene files.'
    usage = 'sgt_constructor.py -i path/to/input/ [OPTIONS]'
    parser, optional, required = help_formatter.initialize_argparse(name='sgt_constructor.py',
                                                                    desc=description,
                                                                    usage=usage)

    # Optional Arguments
    optional.add_argument('-t', '--threads', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Total number of cores/threads for Snakemake to use for parallel job execution.
                          Default: 1
                          """))
    optional.add_argument('--threads-per-job', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Number of threads each individual tool (mafft, raxml-ng) should use.
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
    optional.add_argument('--length_filter_only', action='store_true',
                          help=textwrap.dedent("""\
                          Only perform length filtration.
                          No final alignment, trimming, or tree building.
                          """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Format of the input files.
                          Options: fasta, phylip (names truncated at 10 characters),
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))
    optional.add_argument('--profile', metavar='<profile>', type=str, default=None,
                          help=textwrap.dedent("""\
                          Snakemake cluster profile to use for running jobs.
                          Default: None
                          """))
    optional.add_argument('-n', '--dry-run', action='store_true',
                          help=textwrap.dedent("""\
                          Perform a dry-run of the Snakemake workflow without executing any jobs.
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

    run_snakemake()
