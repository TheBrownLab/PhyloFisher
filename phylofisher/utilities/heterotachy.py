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
from ete3 import Tree

from phylofisher import help_formatter


def bash(cmd):
    """
    Runs bash commands as a subproccess
    @param cmd:
    @return: NONE
    """
    subprocess.run(cmd, executable='/bin/bash', shell=True)


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

    SeqIO.write(records, 'renamed.phy', 'phylip')

    return pseudonames, pseudonames_rev


def fake_tree(treefile, pseudonames):
    with open('renamed.tre', 'w') as res:
        original = open(treefile).readline()
        for key, value in pseudonames.items():
            original = original.replace(key, value)
        res.write(original)


def get_taxa():
    """
    Parses input supermatrix to retrieve taxa
    @return: list of taxa and a dictionary with names as keys and as pandas series of sites as the keys
    """
    taxa_list = []
    mat_dict = {}
    with open('renamed.phy', 'r') as infile:
        for record in SeqIO.parse(infile, 'phylip'):
            taxa_list.append(record.description)
            mat_dict[record.name] = pd.Series(list(record.seq))
    return taxa_list, mat_dict


def get_branch_lens():
    """
    Parse input tree and retrieve branch lengths
    @return:
    """
    tree = Tree('renamed.tre')
    dist_mat = []
    for i, taxon1 in enumerate(taxa):
        dist_mat.append([])
        node1 = tree & taxon1
        for taxon2 in taxa:
            node2 = tree & taxon2
            dist = tree.get_distance(node1, node2)
            dist_mat[i].append(dist)
    length = int((len(taxa) - 1) / 2)
    taxa_dict = {}
    for i, taxon in enumerate(taxa):
        bran_lens = sorted(dist_mat[i], reverse=True)[0:length]
        taxa_dict[taxon] = sum(bran_lens) / length
    series = pd.Series(taxa_dict)
    size = int(len(taxa) / 3)

    return series, size


def prune_tree(speed):
    """

    @param speed:
    @return:
    """
    # Prune tree into SLOW tree
    if speed == 'fast':
        mylist = list(series.nlargest(size).index)
    else:
        mylist = list(series.nsmallest(size).index)
    tree = Tree('renamed.tre')
    nodes = [tree & x for x in mylist]
    tree.prune(nodes, preserve_branch_length=True)
    tree.write(outfile=f'{speed}.tre')

    return mylist


def trim_matrix():
    """

    @return:
    """
    with open('renamed.phy', 'r') as infile, open('fast.phy', 'w') as fast_fas, open('slow.phy', 'w') as slow_fas:
        fast_recs = []
        slow_recs = []
        for record in SeqIO.parse(infile, 'phylip'):
            if record.description in slow_taxa:
                slow_recs.append(record)
            elif record.description in fast_taxa:
                fast_recs.append(record)

        SeqIO.write(slow_recs, slow_fas, 'phylip')
        SeqIO.write(fast_recs, fast_fas, 'phylip')


def get_site_rates():
    """

    @return:
    """
    for speed in ['slow', 'fast']:
        with open(f'{speed}.dist_est.ctl', 'w') as outfile:
            ctl = (f'treefile = {speed}.tre * treefile\n'
                   f'seqfile = {speed}.phy * sequence data\n'
                   'nchar = 20             * amino acid data\n'
                   'model =  3             * empirical + F\n'
                   'aaRatefile = lg.dat * JTT substitution model\n'

                   'nrate = 101            * number of rates\n'
                   'ub = 10.0              * upper bound for rates\n')
            outfile.write(ctl)

    # Calculate site rates for FAST taxa
    bash(f'dist_est fast.dist_est.ctl')
    os.rename('rate_est.dat', 'fast.rate_est.dat')
    os.rename('DE.dat', 'fast.DE.dat')

    # Calculate site rates for SLOW taxa
    bash(f'dist_est slow.dist_est.ctl')
    os.rename('rate_est.dat', 'slow.rate_est.dat')
    os.rename('DE.dat', 'slow.DE.dat')


def parse_rates(speed):
    """

    @param speed:
    @return:
    """
    if speed == 'fast':
        rate_est = 'fast.rate_est.dat'
    else:
        rate_est = 'slow.rate_est.dat'

    rates = []
    for line in open(rate_est):
        _, myrate, _, _ = line.split()
        rates.append(myrate)

    return rates


def get_site_ratio():
    """

    @return:
    """
    rates_ratio = [float(x) / float(y) for x, y in zip(fast_rates, slow_rates)]
    positions = []
    for i, rate in enumerate(rates_ratio):
        positions.append((i, float(rate)))
    sorted_positions = sorted(positions, key=lambda x: x[1], reverse=True)
    result = [i[0] for i in sorted_positions]

    return result


def site_removal():
    """

    @return:
    """
    os.mkdir(f'steps_{args.step_size}')
    os.chdir(f'steps_{args.step_size}')

    out_dict = {'fasta': 'fas',
                'phylip': 'phy',
                'phylip-relaxed': 'phy',
                'nexus': 'nex'}

    for i, step in enumerate(range(args.step_size, len(sorted_sites), args.step_size)):
        with open(f'step{i}.{out_dict[args.in_format]}', 'w') as outfile:
            records = []
            for name, seq in matrix_dict.items():
                seq = "".join(seq[sorted_sites[step:]].values)
                records.append(SeqRecord(Seq(seq),
                                         id=name_dict_rev[name],
                                         name='',
                                         description=''))

            # Writes to output matrix in user specified output
            if args.out_format.lower() in out_dict:
                SeqIO.write(records, outfile, args.out_format.lower())
            else:
                sys.exit('Invalid Output Format')


if __name__ == '__main__':
    description = ('Removes the most heterotachous sites within a phylogenomic supermatrix in a stepwise fashion,'
                   ' leading to a user defined set of new matrices with these sites removed.')
    usage = 'heterotachy.py -t path/to/tree -m path/to/matrix [OPTIONS]'
    parser, optional, required = help_formatter.initialize_argparse(name='heterotachy.py',
                                                                    desc=description,
                                                                    usage=usage)

    required.add_argument('-tr', '--tree', type=str, metavar='', required=True,
                          help=textwrap.dedent("""\
                          Path to tree
                          """))
    required.add_argument('-m', '--matrix', type=str, metavar='', required=True,
                          help=textwrap.dedent("""\
                          Path to supermatrix
                          """))
    # Optional Arguments

    optional.add_argument('-s', '--step_size', type=int, default=3000, metavar='N',
                          help=textwrap.dedent("""\
                          Size of removal step (i.e., 1000 sites removed) to exhaustion
                          Default: 3000
                          """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Input format of matrix
                          Options: fasta, nexus, phylip (names truncated at 10 characters), 
                          or phylip-relaxed (names are not truncated)
                          Default: phylip-relaxed
                          """))
    optional.add_argument('-f', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Desired format of the output steps.
                          Options: fasta, nexus, phylip (names truncated at 10 characters), 
                          or phylip-relaxed (names are not truncated)
                          Default: phylip-relaxed
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

    cwd = os.getcwd()

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    os.chdir(args.output)

    name_dict, name_dict_rev = fake_phylip(f'{cwd}/{args.matrix}')
    fake_tree(f'{cwd}/{args.tree}', name_dict)

    # Parse input alignment
    taxa, matrix_dict = get_taxa()

    # Get Branch Lengths
    series, size = get_branch_lens()

    # Prune tree
    slow_taxa = prune_tree('slow')
    fast_taxa = prune_tree('fast')

    # Slow and Fast fastas
    trim_matrix()

    get_site_rates()
    fast_rates = parse_rates('fast')
    slow_rates = parse_rates('slow')

    sorted_sites = get_site_ratio()
    site_removal()
