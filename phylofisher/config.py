#!/usr/bin/env python
import argparse
import configparser
import textwrap
from pathlib import Path

from phylofisher import help_formatter

if __name__ == '__main__':
    description = 'Script for the analysis folder configuration.'
    usage = 'config.py -d <database_folder> -i <input_file.tsv> [OPTIONS]'
    parser, optional, required = help_formatter.initialize_argparse(name='config.py',
                                                                    desc=description,
                                                                    usage=usage)

    # Required Arguments
    required.add_argument('-d', '--database_folder', required=True, metavar='<database_dir>',
                          help=textwrap.dedent("""\
                          Path to directory containing the PhyloFisher database.
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
                              Path to alternative single gene tree color configuration file.
                              """))

    config = configparser.ConfigParser()
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    if not args.orthomcl:
        args.orthomcl = str(Path(args.database_folder, 'orthomcl'))

    if not args.tree_colors:
        args.tree_colors = str(Path(args.database_folder, 'tree_colors.csv'))

    with open('config.ini', 'w') as configfile:
        config['PATHS'] = {'database_folder': args.database_folder,
                           'input_file':     args.input_file,
                           'orthomcl':       args.orthomcl,
                           'color_conf':     args.tree_colors}
        config.write(configfile)
