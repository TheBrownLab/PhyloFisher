#!/usr/bin/env python
import os
from Bio import SeqIO
import glob
import argparse

def good_length(trimmed_aln):
    records = list(SeqIO.parse(trimmed_aln, 'fasta'))
    with open(trimmed_aln, 'w') as res:
        for record in records:
            coverage = len(str(record.seq).replace('-', '').replace('X', '')) / len(record.seq)
            if coverage > 0.3:
                res.write(f'>{record.description}_{round(coverage, 2)}\n{record.seq}\n')
            else:
                print(record.name, trimmed_aln, coverage)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scriptm which filters short sequences and adds'
                                                 'information about lenght to all other seqs.',
                                     usage="len_filter -i <folder_with_trimmed_alignments> [options]")
    parser.add_argument('-i', '--input', required=True, help='folder with alignments')
    parser.add_argument('-s', '--sufix', help='Sufix of your alignments. Not required is folder'
                                             'contains only alignments.')
    args = parser.parse_args()

    if args.sufix:
        files = glob.glob(os.path.join(args.input,f'*{args.sufix}'))
    else:
        files = glob.glob(os.path.join(args.input,'*'))

    for file in files:
        good_length(file)
