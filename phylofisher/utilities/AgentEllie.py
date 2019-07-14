#!/usr/bin/env python
import os
import subprocess
import argparse
from Bio import SeqIO
import pandas as pd
import string
import random


def id_generator(size=10, chars=string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def unique_name(keys):
    id_ = id_generator()
    if id_ not in keys:
        return id_
    else:
        unique_name(keys)

def fake_phylip(matrix):
    seqs = 0
    length = None
    pseudonames = {}
    pseudonames_rev = {}
    result = ''
    with open('TEMP.phy', 'w') as res:
        for record in SeqIO.parse(matrix, 'fasta'):
            seqs += 1
            if not length:
                length = len(record.seq)
            uname = unique_name(pseudonames)
            pseudonames[record.name] = uname
            pseudonames_rev[uname] = record.name
            result += f'{uname} {record.seq}\n'
        res.write(f' {seqs} {length}\n')
        res.write(result)
    return pseudonames, pseudonames_rev

def fake_tree(treefile, pseudonames):
    with open('TEMP.tre', 'w') as res:
        original = open(treefile).readline()
        for key, value in pseudonames.items():
            original = original.replace(key, value)
        res.write(original)

def control_file():
    ctl = """treefile = TEMP.tre * treefile"
seqfile = TEMP.phy * sequence data

nchar = 20             * amino acid data
model =  3             * empirical + F
aaRatefile = lg.dat * JTT substitution model

nrate = 101            * number of rates
ub = 10.0              * upper bound for rates"""
    with open('dist_est.ctl', 'w') as res:
        res.write(ctl)


def run_dist():
    cmd = './dist_est dist_est.ctl'
    subprocess.run(cmd, shell=True)


def parse_rates():
    positions = []
    site = 0
    for line in open('rate_est.dat'):
        _, rate, _, _ = line.split()
        positions.append((site, float(rate)))
        site += 1
    sorted_positions = sorted(positions, key=lambda x: x[1], reverse=True)
    result = [i[0] for i in sorted_positions]
    return result


def main():
    pseudo_, pseudo_rev_ = fake_phylip(args.matrix)
    fake_tree(args.tree, pseudo_)
    control_file()
    run_dist()
    matrix_dict = {}
    for record in SeqIO.parse(args.matrix, 'fasta'):
        matrix_dict[record.name] = pd.Series(list(record.seq))

    sorted_rates = parse_rates()
    iter = 0
    os.mkdir(f'chunks_{args.chunk}')
    os.chdir(f'chunks_{args.chunk}')
    for chunk in range(args.chunk, len(sorted_rates), args.chunk):
        with open(f'chunk{iter}', 'w') as res:
            res.write(f' {len(matrix_dict)} {len(sorted_rates) - chunk}\n')
            for name, seq in matrix_dict.items():
                res.write(f'{name} {"".join(seq[chunk:].values)}\n')
        iter += 1
    os.remove('../TEMP.phy')
    os.remove('../TEMP.tre')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='some description', usage="blabla")
    parser.add_argument('-m', '--matrix')
    parser.add_argument('-tr', '--tree')
    parser.add_argument('-c', '--chunk', type=int, default=3000)
    args = parser.parse_args()
    main()









