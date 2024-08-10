#!/usr/bin/env python
import os
import random
import string
import subprocess
import sys
import textwrap

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofisher import help_formatter

# Modified from https://github.com/TheBrownLab/PhyloFisher/blob/master/phylofisher/utilities/fast_site_remover.py

def id_generator(size=10, chars=string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def unique_name(keys):
    id_ = id_generator()
    if id_ not in keys:
        return id_
    else:
        unique_name(keys)


def fake_phylip(matrix):
    length = None
    pseudonames = {}
    pseudonames_rev = {}
    records = []

    for record in SeqIO.parse(matrix, args.in_format):
        uname = unique_name(pseudonames)
        pseudonames[record.name] = uname
        pseudonames_rev[uname] = record.name
        seq = str(record.seq)
        records.append(SeqRecord(Seq(seq),
                                 id=uname,
                                 name='',
                                 description=''))

    SeqIO.write(records, 'TEMP.phy', 'phylip')

    return pseudonames, pseudonames_rev


def fake_tree(treefile, pseudonames):
    with open('TEMP.tre', 'w') as res:
        original = open(treefile).readline()
        for key, value in pseudonames.items():
            original = original.replace(key, value)
        res.write(original)

def control_file_nuc():
    ctl = """treefile = TEMP.tre * treefile
    seqfile = TEMP.phy * sequence data
    
    nchar = 4           * nucleic data
    model = 2           * F81 (don't want to give rates for HKY/GTR...)
    nrate = 101         * number of rates
    ub = 10             * upper bound for rates"""
    with open('dist_est.ctl', 'w') as res:
        res.write(ctl)

def control_file_amino():
    ctl = """treefile = TEMP.tre * treefile
seqfile = TEMP.phy * sequence data

nchar = 20             * amino acid data
model =  3             * empirical + F
aaRatefile = lg.dat * JTT substitution model

nrate = 101            * number of rates
ub = 10.0              * upper bound for rates"""
    with open('dist_est.ctl', 'w') as res:
        res.write(ctl)


def run_dist():
    cmd = 'dist_est dist_est.ctl'
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
    if args.type == "amino":
        control_file_amino()
    elif args.type == "nuc":
        control_file_nuc()
    else:
        ValueError("Type should be amino or nuc")
    # Checks to see if dist_est has been run. If not do so
    if not os.path.isfile('./rate_est.dat'):
        run_dist()

    matrix_dict = {}
    for record in SeqIO.parse(args.matrix, args.in_format):
        matrix_dict[record.name] = pd.Series(list(record.seq))

    sorted_rates = parse_rates()
    iter = 0
    os.mkdir(f'{args.output}/steps_{args.step_size}')
    os.chdir(f'{args.output}/steps_{args.step_size}')

    if args.pos == True:
        with open(f'fastest_sites.txt', 'w') as res:
            res.write(str(sorted_rates))

    out_dict = {'fasta'         : 'fas',
                'phylip'        : 'phy',
                'phylip-relaxed': 'phy',
                'nexus'         : 'nex'}

    for step in range(args.step_size, len(sorted_rates), args.step_size):
        with open(f'step{iter}.{out_dict[args.out_format.lower()]}', 'w') as res:
            records = []
            for name, seq in matrix_dict.items():
                seq = "".join(seq[sorted_rates[step:]].values)
                records.append(SeqRecord(Seq(seq),
                                         id=name,
                                         name='',
                                         description=''))

            # Writes to output matrix in user specified output
            if args.out_format.lower() in out_dict:
                SeqIO.write(records, res, args.out_format.lower())
            else:
                sys.exit('Invalid Output Format')

        iter += 1
    os.remove('../../TEMP.phy')
    os.remove('../../TEMP.tre')


if __name__ == "__main__":
    description = 'This removes the fastest evolving sites within the phylogenomic supermatrix in a stepwise fashion.'
    parser, optional, required = help_formatter.initialize_argparse(name='fast_site_remover.py',
                                                                    desc=description,
                                                                    usage='fast_site_remover.py '
                                                                          '[OPTIONS] -i /path/to/input/')

    # Required Arguments
    required.add_argument('-m', '--matrix', required=True, type=str, metavar='<matrix>',
                          help=textwrap.dedent("""\
                              Path to matrix
                              """))
    required.add_argument('-tr', '--tree', type=str, metavar='',
                          help=textwrap.dedent("""\
                                Path to tree
                                  """))

    # Optional Arguments
    optional.add_argument('-s', '--step_size', type=int, default=3000, metavar='N',
                          help=textwrap.dedent("""\
                                Size of removal step (i.e., 1000 sites removed) to exhaustion
                                Default: 3000
                                """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                                Format of input matrix.
                                Options: fasta, nexus, phylip (names truncated at 10 characters), 
                                or phylip-relaxed (names are not truncated)
                                Default: fasta
                                """))
    optional.add_argument('-of', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                                Desired format of the output steps.
                                Options: fasta, nexus, phylip (names truncated at 10 characters), 
                                or phylip-relaxed (names are not truncated)
                                Default: fasta
                                """))
    optional.add_argument('-t', '--type', type=str, default='amino',
                          help=textwrap.dedent("""\
                                               Whether file is amino or nuc (nucleic)
                                               Default: amino
                                               """))
    optional.add_argument('-p', '--pos', type=bool, default=False,
                          help=textwrap.dedent("""\
                                Prints out fastest site positions to a file fastest_sites.txt in order of fastest
                                               to slowest site if True.
                                Default: False
                                               """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

    os.mkdir(args.output)

    main()
