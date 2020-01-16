#!/usr/bin/env python
import os
import textwrap
from glob import glob
from Bio import SeqIO
from collections import defaultdict
import argparse
import csv


def parse_names(input_folder):
    """

    :param input_folder:
    :return:
    """
    name_set = set()
    if args.suffix:
        files = sorted(glob(f'{input_folder}/*{args.suffix}'))
    else:
        files = sorted(glob(f'{input_folder}/*'))
    for file in files:
        with open(file) as f:
            for line in f:
                if line.startswith('>'):
                    fname = line.split()[0][1:]
                    name = fname.split('_')[0]
                    name_set.add(name)
    return files, sorted(list(name_set))


def stats(total_len):
    """

    :param total_len:
    :return:
    """
    with open(f'{args.output}_forge_stats.tsv', 'w') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['org', 'missing[%]'])
        missing = []
        for record in SeqIO.parse(args.output, 'fasta'):
            missing.append((record.name, (record.seq.count('-')/total_len)*100))
        for org_missing in sorted(missing, key=lambda x: x[1], reverse=True):
            tsv_writer.writerow(list(org_missing))


def main():
    """

    :return: NONE
    """
    input_folder = os.path.basename((args.input.strip('/')))
    files, orgs = parse_names(input_folder)
    total_len = 0
    res_dict = defaultdict(str)
    with open(f'{args.output}_indices.tsv', 'w') as outfile:
        for file in files:
            gene = os.path.basename(file).split('.')[0]
            length = 0
            seq_dict = {}
            for record in SeqIO.parse(file, 'fasta'):
                length = len(record.seq)
                seq_dict[record.id.split('_')[0]] = str(record.seq)
            start_len = total_len + 1
            total_len += length
            outfile.write(f'{gene}\t{start_len}\t{total_len}\n')
            for org in orgs:
                if org in seq_dict:
                    res_dict[org] += seq_dict[org]
                else:
                    res_dict[org] += ('-' * length)

        with open(args.output, "w") as res:
            # Writes in different output formats

            # fasta format
            if args.format.lower() is 'fasta':
                for org, seq in res_dict.items():
                    res.write(f'>{org}\n{seq}\n')

            # phylip format TODO: Test PHYLIP formatted ouput
            elif args.format.lower() is 'phylip':
                for org, seq in res_dict.items():
                    # Makes seq id == exactly 10 characters
                    if len(org) <= 10:
                        org = org + ((10 - len(org)) * ' ')
                    else:
                        org = org[0:9]
                    res.write(f'{org}{seq}\n')

            # nexus format TODO: Test NEXUS formatted output
            elif args.format.lower() is 'nexus':
                entries = ''
                count = 0
                for org, seq in res_dict.items():
                    entries += f'      {org}\t{seq}\n'
                    count += 1
                res.write(f'#nexus\n'
                          f'...\n'
                          f'begin data;\n'
                          f'dimensions ntax={count} nchar={total_len};\n'
                          f'   format datatype=dna missing=? gap=-;\n'
                          f'   matrix\n'
                          f'{entries}'
                          f'   ;\n'
                          f'end;\n')

            # TODO: Raise an exception for "invalid output type"
            else:
                pass
        stats(total_len)


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
                                     description='some description',
                                     usage='forge.py [OPTIONS] -i /path/to/input/',
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
    optional.add_argument('-o', '--output', default="output", type=str, metavar='',
                          help=textwrap.dedent("""\
                          Desired basename of output files. 
                          Default: output
                          Example: output, output_indices.tsv, & output_forge_stats.tsv
                          """))
    optional.add_argument('-f', '--out_format', metavar='"format"', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Desired format of the output matrix.
                          Options: fasta, phylip, or nexus.
                          Default: fasta
                          """))
    optional.add_argument('-s', '--suffix', metavar='"suffix"', type=str,
                          help=textwrap.dedent("""\
                          Suffix of input files
                          Default: NONE
                          Example: path/to/input/*.suffix
                          """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                          Show this help message and exit.
                          """))

    parser._action_groups.append(optional)
    args = parser.parse_args()

    main()