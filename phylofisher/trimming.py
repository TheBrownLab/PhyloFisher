#!/usr/bin/env python
import os
import argparse
import textwrap
from glob import glob
from pathlib import Path


# TODO fucking cleaning


def prepare_analyses(dataset, threads):
    root = dataset.split('/')[-1].split('.')[0]
    command = (f'no_gap_stops.py {dataset} &&'
               
               f'prequal {root}.aa && '
               
               f'mafft --globalpair --maxiterate 1000 --unalignlevel 0.6'
               f' --thread {threads} {root}.aa.filtered > {root}.aln && '
               
               f'divvier -mincol 4 -divvygap {root}.aln && '
               
               f'pre_trimal.py {root}.aln.divvy.fas && '
               
               f'BMGE -t AA -g 0.3 -i {root}.pre_trimal -of {root}.bmge && '
               
               f'len_filter2.py -i {root}.bmge -t 0.5 && '
               
               f'mafft --globalpair --maxiterate 1000 --unalignlevel 0.6'
               f' --thread {threads} {root}.len > {root}.aln2 && '
               
               f'divvier -mincol 4 -divvygap {root}.aln2 && '
               
               f'trimal -in {root}.aln2.divvy.fas -gt 0.01 -out {root}.final\n'
               
               f'raxmlHPC-PTHREADS-AVX2 -T {threads} -m PROTGAMMALG4XF -f a -s {root}.final'
               f' -n {root}.tre -x 123 -N 100 -p 12345 && '
               
               f'add_aln_length.py {root}\n')

    return command


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
    parser = argparse.ArgumentParser(prog='trimming.py',
                                     description='description',
                                     usage='trimming.py -i path/to/input/ [OPTIONS]',
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

    # Optional Arguments
    optional.add_argument('-s', '--suffix', metavar='"suffix"', type=str,
                          help=textwrap.dedent("""\
                              Suffix of input files
                              Default: NONE
                              Example: path/to/input/*.suffix
                              """))
    optional.add_argument('-t', '--threads', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Desired number of threads to be utilized.
                          Default: 1
                          """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                          Show this help message and exit.
                          """))

    parser._action_groups.append(optional)
    args = parser.parse_args()

    lib = f'{os.path.realpath(__file__).split("phylofisher")[0]}lib'

    os.chdir(args.input_folder)
    with open('commands.txt', 'w') as res:
        for file in glob(f'*{args.suffix}'):
            line = prepare_analyses(file, args.threads)
            res.write(line)
