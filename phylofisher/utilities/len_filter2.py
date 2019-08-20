#!/usr/bin/env python
import os
from Bio import SeqIO
import argparse

def read_full_proteins(core):
    full_prots = {}
    for record in SeqIO.parse(f'{core}.aa', 'fasta'):
        full_prots[record.name] = record.seq
    return full_prots


def good_length(trimmed_aln, threshold):
    core = trimmed_aln.split('.')[0]
    full_proteins = read_full_proteins(core)
    original_name = f'{core}.len'
    length = None
    with open(original_name, 'w') as res:
        for record in SeqIO.parse(trimmed_aln, 'fasta'):
            if length is None:
                length = len(record.seq)
            coverage = len(str(record.seq).replace('-', '').replace('X', '')) / len(record.seq)
            if coverage > threshold:
                res.write(f'>{record.description}_{round(coverage, 2)}\n{full_proteins[record.name]}\n')
            else:
                print('deleted:', record.name, coverage)
    with open(f'{core}.length', 'w') as f:
        f.write(length)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script which filters short sequences and adds'
                                                 'information about lenght to all other seqs.',
                                     usage="len_filter -i <input_alignment> -t <len threshold>")
    parser.add_argument('-i', '--input', required=True, help='trimmed alignment input')
    parser.add_argument('-t', '--threshold', default=0.5, type=float, help='len threshold e.g. 0.5 for 50%%')
    args = parser.parse_args()
    good_length(args.input, args.threshold)
