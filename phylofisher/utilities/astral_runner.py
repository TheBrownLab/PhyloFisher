import os
import shutil
import subprocess
import textwrap

from phylofisher import help_formatter


def make_astral_inputs():
    """

    :return:
    """
    genes = [file.split('.')[1] for file in os.listdir(args.input) if 'RAxML_bipartitions.' in file]
    genes = set(genes)
    
    with open(f'{args.output}/sgt_trees.tre', 'w') as tree_file, open(f'{args.output}/bs-files.txt', 'w') as bs_files:
        for gene in genes:
            with open(f'{args.input}/RAxML_bipartitions.{gene}.{args.suffix}', "r") as infile:
                for line in infile:
                    line = line.strip()
                    tree_file.write(f'{line}\n')
            
            bs_files.write(f'{args.input}/RAxML_bootstrap.{gene}.{args.suffix}\n')
                    

def run_astral():
    """

    :return:
    """
    os.chdir(args.output)
    subprocess.run(f'astral -i sgt_trees.tre -b bs-files.txt -o AstralBS.out', shell=True, executable='/bin/bash')


if __name__ == '__main__':
    description = 'Description'  # TODO: Get Description
    parser, optional, required = help_formatter.initialize_argparse(name='astral_runner.py',
                                                                    desc=description,
                                                                    usage='astral_runner.py '
                                                                          '[OPTIONS] -i /path/to/input/')

    # Optional Arguments
    # optional.add_argument('--in_format', metavar='<format>', type=str, default='phylip-relaxed',
    #                       help=textwrap.dedent("""\
    #                               Desired format of the output chunks.
    #                               Options: fasta, nexus, phylip (names truncated at 10 characters),
    #                               or phylip-relaxed (names are not truncated)
    #                               Default: phylip-relaxed
    #                               """))

    in_help = 'Path to input directory containing gene files in FASTA format.'  # TODO: Correct this
    args = help_formatter.get_args(parser, optional, required, in_help=in_help)
    
    args.output = os.path.abspath(args.output)
    args.input = os.path.abspath(args.input)
    
    if os.path.isdir(args.output):
        shutil.rmtree(args.output)
    os.mkdir(args.output)

    make_astral_inputs()
    run_astral()
