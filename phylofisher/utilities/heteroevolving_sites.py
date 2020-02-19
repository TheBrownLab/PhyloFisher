#!/usr/bin/env python
import os
from ete3 import Tree
from Bio import SeqIO
import argparse
from phylofisher.utilities.fast_tax_removal import Leaves


def fast_and_slow(matrix, format, portion, sorted_orgs):
    amount = round(len(sorted_orgs)*portion)
    slow = set(sorted_orgs[(len(sorted_orgs)-amount):])
    fast = set(sorted_orgs[:amount])
    print(len(sorted_orgs))
    print(len(slow))
    print(len(fast))
    with open('slow.fas', 'w') as s, open('fast.fas', 'w') as f:
        for record in SeqIO.parse(matrix, format):
            if record.name in slow:
                s.write(f'>{record.name}\n{record.seq}\n')
            elif record.name in fast:
                f.write(f'>{record.name}\n{record.seq}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fast taxon removal', usage="fast_tax_removal.py [OPTIONS]")
    parser.add_argument('-t', '--tree', required=True)
    parser.add_argument('-m', '--matrix', required=True)
    parser.add_argument('-f', '--format', default='fasta', help='format of your matrix [default: fasta]')
    parser.add_argument('-p', '--portion', default=0.2, type=float)
    args = parser.parse_args()

    tree = Tree(args.tree)
    leaves = Leaves(tree)
    sorted_taxa = leaves.org_speed()
    fast_and_slow(args.matrix, args.format, args.portion, sorted_taxa)