#!/usr/bin/env python
from Bio import SeqIO
import argparse


def good_length(trimmed_aln, threshold):
    with open(trimmed_aln + '.len', 'w') as res:
        for record in SeqIO.parse(trimmed_aln, 'fasta'):
            coverage = len(str(record.seq).replace('-', '').replace('X', '')) / len(record.seq)
            if coverage > threshold:
                res.write(f'>{record.description}_{round(coverage, 2)}\n{record.seq}\n')
            else:
                print('deleted:', record.name, coverage)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scriptm which filters short sequences and adds'
                                                 'information about lenght to all other seqs.',
                                     usage="len_filter -i <input_alignment> -t <len threshold>")
    parser.add_argument('-i', '--input', required=True, help='  trimmed alignment input')
    parser.add_argument('-t', '--threshold', required=True, type=float, help='len threshold e.g. 0.3 for 30%')
    args = parser.parse_args()
    good_length(args.input, args.threshold)
