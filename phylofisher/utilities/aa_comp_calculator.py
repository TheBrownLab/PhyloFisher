#!/usr/bin/env python
import argparse
import textwrap

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc
from phylofisher import help_formatter


def make_plot():
    df = pd.read_csv(f'{args.output}.tsv', sep="\t")
    df = df.set_index('Taxon')

    # Todo: Make plot prettier
    plt.figure(figsize=(25, 10))
    plt.title("AA Composition Hierarchical Clustering")
    plt.xticks(rotation=90)
    shc.dendrogram(shc.linkage(df, method='ward'), labels=df.index.values, color_threshold=0)
    plt.tight_layout()
    plt.savefig(f'AA_Composition_Hierarchical_Clustering.pdf')


if __name__ == '__main__':
    description = 'Calculates amino acid composition of the supplied super matrix'
    parser, optional, required = help_formatter.initialize_argparse(name='aa_comp_calculator.py',
                                                                    desc=description,
                                                                    usage='aa_comp_calculator.py '
                                                                          '[OPTIONS] -i <matrix>')

    required.add_argument('-i', '--input', required=True, type=str, metavar='matrix',
                          help=textwrap.dedent("""\
                                  Path to input matrix for analysis.
                                  """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

    peptides = ['A', 'G', 'P', 'S', 'T', 'C', 'F', 'W', 'Y', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'N', 'D', 'E', 'Q']

    with open(args.input, 'r') as infile, open(f'{args.output}.tsv', 'w') as outfile:
        outfile.write('Taxon\t'+'\t'.join(peptides) + '\n')

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
