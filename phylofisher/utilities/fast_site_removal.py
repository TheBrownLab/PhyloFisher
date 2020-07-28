#!/usr/bin/env python
import os
import random
import string
import subprocess
import sys
import textwrap

import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofisher import help_formatter


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


def write_seqs(out_handle, records):
    # Accepted out formats with respective suffix
    out_dict = {'fasta': 'fas',
                'phylip': 'phy',
                'phylip-relaxed': 'phy',
                'nexus': 'nex'}

    # Writes to output matrix in user specified output
    if args.out_format.lower() in out_dict:
        SeqIO.write(records, out_handle, args.out_format.lower())
    else:
        sys.exit('Invalid Output Format')


def main():
    pseudo_, pseudo_rev_ = fake_phylip(args.matrix)
    fake_tree(args.tree, pseudo_)
    control_file()
    # Checks to see if dist_est has been run. If not do so
    if not os.path.isfile('./rate_est.dat'):
        run_dist()

    matrix_dict = {}
    for record in SeqIO.parse(args.matrix, 'fasta'):
        matrix_dict[record.name] = pd.Series(list(record.seq))

    sorted_rates = parse_rates()
    iter = 0
    os.mkdir(f'{args.output}/chunks_{args.chunk}')
    os.chdir(f'{args.output}/chunks_{args.chunk}')
    for chunk in range(args.chunk, len(sorted_rates), args.chunk):
        with open(f'chunk{iter}', 'w') as res:
            records = []
            for name, seq in matrix_dict.items():
                seq = "".join(seq[sorted_rates[chunk:]].values)
                records.append(SeqRecord(Seq(seq, IUPAC.protein),
                                         id=name,
                                         name='',
                                         description=''))
            write_seqs(res, records)

        iter += 1
    # os.remove('../TEMP.phy')
    # os.remove('../TEMP.tre')


if __name__ == "__main__":
    description = 'This removes the fastest evolving sites within the phylogenomic supermatrix in a stepwise fashion.'
    parser, optional, required = help_formatter.initialize_argparse(name='fast_site_removal.py',
                                                                    desc=description,
                                                                    usage='fast_site_removal.py '
                                                                          '[OPTIONS] -i /path/to/input/')

    # Required Arguments
    required.add_argument('-m', '--matrix', required=True, type=str, metavar='<matrix>',
                          help=textwrap.dedent("""\
                              Path to matrix
                              """))

    # Optional Arguments
    optional.add_argument('-tr', '--tree', type=str, metavar='',
                          help=textwrap.dedent("""\
                              Path to tree
                              """))
    optional.add_argument('-c', '--chunk', type=int, default=3000, metavar='N',
                          help=textwrap.dedent("""\
                              Size of removal step (i.e., 1000 sites removed) to exhaustion
                              Default: 3000
                              """))
    optional.add_argument('-f', '--out_format', metavar='<format>', type=str, default='phylip-relaxed',
                          help=textwrap.dedent("""\
                              Desired format of the output chunks.
                              Options: fasta, nexus, phylip (names truncated at 10 characters), 
                              or phylip-relaxed (names are not truncated)
                              Default: phylip-relaxed
                              """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)
    
    os.mkdir(args.output)
    
    main()
