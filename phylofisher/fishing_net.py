#!/usr/bin/env python
import os
import argparse
import textwrap
from glob import glob
import configparser
from pathlib import Path
from Bio import SeqIO


def parse_genes(gene_file):
    to_exlude = set()
    with open(gene_file) as lines:
        next(lines)
        for line in lines:
            gene, _, _, sgt = line.split(',')
            sgt = sgt.strip()
            if sgt.lower() != 'yes':
                to_exlude.add(gene)
    return to_exlude


def parse_orgs(org_file):
    to_exlude = set()
    paralogs = set()
    with open(org_file) as lines:
        header = next(lines)
        if header.count(',') == 10:
            sgt_idx = 9
        else:
            sgt_idx = 6
        for line in lines:
            org = line.split(',')[0]
            sgt = line.split(',')[sgt_idx].strip()
            if sgt.lower() != 'yes':
                to_exlude.add(org)
            if header.count(',') == 10:
                if line.split(',')[sgt_idx + 1].strip().lower() == "yes":
                    if org not in to_exlude:
                        paralogs.add(org)
    return to_exlude, paralogs


def fasta_filtr(file, o_to_ex, paralogs=None):
    with open(str(Path(args.output_directory, file)), 'w') as res:
        for record in SeqIO.parse(str(Path(args.input_directory, file)), 'fasta'):
            if record.name.split('_')[0] not in o_to_ex:
                res.write(f'>{record.name}\n{record.seq}\n')
        if paralogs:
            para_file = str(Path(dfo, f'paralogs/{file.split(".")[0]}_paralogs.fas'))
            if os.path.isfile(para_file):
                for record in SeqIO.parse(para_file, 'fasta'):
                    if record.name.split('.')[0] in paralogs:
                        res.write(f'>{record.name}\n{record.seq}\n')


def main():
    # if args.orthologs:
    #     gene_file = 'orthologs_stats/genes_stats.csv'
    #     orgs_file = 'orthologs_stats/orgs_stats.csv'
    # else:
    gene_file = str(Path(args.input_directory + '_stats', 'genes_stats.csv'))
    orgs_file = str(Path(args.input_directory + '_stats', 'orgs_stats.csv'))
    g_to_ex = parse_genes(gene_file)
    o_to_ex, paralogs = parse_orgs(orgs_file)
    filtered_genes = []
    for file in glob(args.input_directory + f'/*{args.sufix}'):
        if file.split('.')[0] not in g_to_ex:
            filtered_genes.append(file)

    for file in filtered_genes:
        if paralogs:
            fasta_filtr(os.path.basename(file), o_to_ex, paralogs)
        else:
            fasta_filtr(os.path.basename(file), o_to_ex)


if __name__ == '__main__':
    class CustomHelpFormatter(argparse.HelpFormatter):
        """This class can be used to make visual changes in the help"""

        def _format_action_invocation(self, action):
            # This removes metvar after short option
            if not action.option_strings or action.nargs == 0:
                return super()._format_action_invocation(action)
            default = self._get_default_metavar_for_optional(action)
            args_string = self._format_args(action, default)
            return ', '.join(action.option_strings) + ' ' + args_string

        def _split_lines(self, text, width):
            # This adds 3 spaces before lines that wrap
            lines = text.splitlines()
            for i in range(0, len(lines)):
                if i >= 1:
                    lines[i] = (3 * ' ') + lines[i]
            return lines


    class myHelpFormatter(CustomHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    formatter = lambda prog: myHelpFormatter(prog, max_help_position=100)
    parser = argparse.ArgumentParser(prog='forge.py',
                                     description='Script for filtering orgs [and|or] genes',
                                     usage='fishing_net.py -i <input_directory> -o <output_directory> [OPTIONS]',
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent("""\
                                     additional information:
                                        stuff
                                        """))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    # Required Arguments
    required.add_argument('-i', '--input', required=True, type=str, metavar='path/to/input/',
                          help=textwrap.dedent("""\
                          Path to input directory
                          """))
    required.add_argument('-o', '--output', default="output", type=str, metavar='',
                          help=textwrap.dedent("""\
                          Path to desired output directory
                          """))

    # Optional Arguments

    optional.add_argument('-d', '--dataset_folder', metavar='<dataset>', type=str, default='fasta',  # TODO: get default
                          help=textwrap.dedent("""\
                          Path to directory containing the dataset
                          """))
    optional.add_argument('-s', '--suffix', metavar='<suffix>', type=str,
                          help=textwrap.dedent("""\
                          Suffix of input files
                          """))
    optional.add_argument('--orthologs', action='store_true',
                          help=textwrap.dedent("""\
                          Only for ortholog selection. Without information
                          about used path."""))
    optional.add_argument('-c', '--use_config', action='store_true',
                          help=textwrap.dedent("""\
                          Use config file
                          """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                          Show this help message and exit.
                          """))

    parser._action_groups.append(optional)
    args = parser.parse_args()

    if args.use_config:
        config = configparser.ConfigParser()
        config.read('config.ini')
        dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
        if args.orthologs:
            args.input_directory = os.path.join(dfo, 'orthologs')

    if args.dataset_folder:
        dfo = str(Path(args.dataset_folder).resolve())

    if args.input_directory[-1] == '/':
        args.input_directory = args.input_directory[:-1]
    os.mkdir(args.output_directory)
    main()
