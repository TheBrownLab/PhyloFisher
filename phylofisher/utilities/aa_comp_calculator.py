#!/usr/bin/env python
import argparse
import textwrap

import pandas as pd
import phylofisher.help_formatter
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc


def get_args():
    formatter = lambda prog: phylofisher.help_formatter.myHelpFormatter(prog, max_help_position=100)

    parser = argparse.ArgumentParser(prog='aa_comp_calculator.py',
                                     description='some description',
                                     usage='aa_comp_calculator.py [OPTIONS] -i /path/to/input/',
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
                              Path to input directory containing trimmed alignments in FASTA format.
                              """))

    # Optional Arguments
    optional.add_argument('-o', '--output', default="output", type=str, metavar='',
                          help=textwrap.dedent("""\
                              Desired name of output tsv files. 
                              Default: output
                              Example: output.tsv hcluster.output.pdf
                              """))
    optional.add_argument('-f', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                              Desired format of the output matrix.
                              Options: fasta, phylip (names truncated at 10 characters), 
                              phylip-relaxed (names are not truncated), or nexus.
                              Default: fasta
                              """))
    optional.add_argument('-plot_df', '--suffix', metavar='<suffix>', type=str,
                          help=textwrap.dedent("""\
                              Suffix of input files
                              Default: NONE
                              Example: path/to/input/*.suffix
                              """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                              Show this help message and exit.
                              """))

    return parser.parse_args()


def make_plot():
    df = pd.read_csv(f'{args.output}.tsv', sep="\t")
    print(df.head(10))

    plt.figure(figsize=(40, 10))
    plt.title("AA Composition Hierarchical Clustering")
    shc.dendrogram(shc.linkage(df, method='ward'), labels=df.index.values)
    plt.tight_layout()
    plt.savefig(f'hcluster.{args.output}.pdf')


if __name__ == '__main__':
    args = get_args()
    peptides = ['A', 'G', 'P', 'S', 'T', 'C', 'F', 'W', 'Y', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'N', 'D', 'E', 'Q']

    with open(args.input, 'r') as infile, open(f'{args.output}.tsv', 'w') as outfile:
        outfile.write('\t'.join(peptides) + '\n')

        # Reads in input file
        for record in SeqIO.parse(infile, format='fasta'):
            outfile.write(f'{record.id}\t')
            analysed_seq = ProteinAnalysis(str(record.seq))
            count_dict = analysed_seq.count_amino_acids()
            length = len(str(record.seq).replace("-", "").replace("X", "").replace("*", ""))
            out_str = ''

            # Loops through peptides and checks to see if it is in count_dict
            for pep in peptides:
                if pep in count_dict.keys():
                    out_str += f'{float(count_dict[pep]) / length}\t'
                else:
                    out_str += '0\t'

            outfile.write(out_str.strip() + '\n')

    make_plot()
