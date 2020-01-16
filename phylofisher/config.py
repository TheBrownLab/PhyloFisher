#!/usr/bin/env python
import configparser
import argparse
import textwrap
from pathlib import Path


class CustomHelpFormatter(argparse.HelpFormatter):
    """This class can be used to make chanes in the help"""

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
                      Path to input file in tsv format
                      """))
# Optional Arguments
optional.add_argument('--orthomcl', metavar='<omcl_data>',
                      help=textwrap.dedent("""\
                      Path to orthomcl if not in dataset_dir
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

with open('config.ini', 'w') as configfile:
    config['PATHS'] = {'dataset_folder': args.dataset_folder,
                       'input_file': args.input_file,
                       'orthomcl': args.orthomcl}
    config.write(configfile)
