import argparse
import textwrap

import phylofisher.help_formatter
from phylofisher import fisher


def get_args():
    global args
    formatter = lambda prog: phylofisher.help_formatter.myHelpFormatter(prog, max_help_position=100)
    parser = argparse.ArgumentParser(prog='Recode_Phylip_SR4classes.py',
                                     # TODO: Get description and usage
                                     description='some description',
                                     usage='some usage',
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent("""\
                                             additional information:
                                             """))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    # TODO: What is optional and required?
    # Required Arguments
    required.add_argument('-i', '--input', required=True, type=str, metavar='',
                          help=textwrap.dedent("""\
                          Input file
                          """))
    required.add_argument('-g', '--output', type=str, required=True, metavar='',
                          help=textwrap.dedent("""\
                          groups
                          """))
    # Optional Aruments
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                          Show this help message and exit.
                          """))
    parser._action_groups.append(optional)
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    # AGPST, C, FWY, HRK, MILV, NDEQ === DAYHOFF classes 0-5
    key = {'A': ['A', 'G', 'N', 'P', 'S', 'T'],
           'C': ['C', 'H', 'W', 'Y'],
           'G': ['D', 'E', 'K', 'Q', 'R'],
           'T': ['F', 'I', 'L', 'M', 'V']
           }

    with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
        outfile.write(infile.readline())

        for line in infile:
            name, seq = line.strip().split()
            for nuc in key.keys():
                for x in key[nuc]:
                    seq = seq.replace(x, nuc)
                    seq = seq.replace("X", "-")

                    outfile.write(name + " " + seq + "\n")
