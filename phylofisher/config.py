#!/usr/bin/env python
import configparser
import argparse
import textwrap
from pathlib import Path
import phylofisher.help_formatter


if __name__ == '__main__':
    formatter = lambda prog: phylofisher.help_formatter.MyHelpFormatter(prog, max_help_position=100)
    parser = argparse.ArgumentParser(prog='forest.py',
                                     description='Script for the analysis folder configuration.',
                                     usage="config.py -d <dataset_folder> -i <input_file.tsv> [OPTIONS]",
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent("""\
                                     additional information:
                                        stuff
                                        """))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    # Required Arguments
    required.add_argument('-d', '--dataset_folder', required=True, metavar='<dataset_dir>',
                          help=textwrap.dedent("""\
                          Path to directory containing dataset.
                          """))
    required.add_argument('-i', '--input_file', required=True, metavar='<input.tsv>',
                          help=textwrap.dedent("""\
                          Path to input metadata in tsv format
                          """))
    # Optional Arguments
    optional.add_argument('--orthomcl', metavar='<omcl_data>',
                          help=textwrap.dedent("""\
                          Path to orthomcl if NOT in dataset_dir
                          """))
    optional.add_argument('--tree_colors', metavar='<omcl_data>',
                          help=textwrap.dedent("""\
                              Path to alternative single genetree color configuration file.
                              """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                          Show this help message and exit.
                          """))

    parser._action_groups.append(optional)
    config = configparser.ConfigParser()
    args = parser.parse_args()

    if not args.orthomcl:
        args.orthomcl = str(Path(args.dataset_folder, 'orthomcl'))

    if not args.tree_colors:
        args.tree_colors = str(Path(args.dataset_folder, 'tree_colors.csv'))

    with open('config.ini', 'w') as configfile:
        config['PATHS'] = {'dataset_folder': args.dataset_folder,
                           'input_file': args.input_file,
                           'orthomcl': args.orthomcl}
        config.write(configfile)
