#!/usr/bin/env python

import sys
import textwrap
import subprocess
import os

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
import pandas as pd
from phylofisher import help_formatter


def bash(cmd):
    """
    Runs bash commands as a subproccess
    @param cmd:
    @return: NONE
    """
    subprocess.run(cmd, executable='/bin/bash', shell=True)


def write_seqs(out_handle, records):
    """
    Writes input sequences to a file
    @param out_handle: file handle to write to
    @param records: list of Seq records
    @return: NONE
    """
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


def get_taxa():
    """
    Parses input supermatrix to retrieve taxa
    @return: list of taxa and a dictionary with names as keys and as pandas series of sites as the keys
    """
    taxa_list = []
    mat_dict = {}
    with open(args.matrix, 'r') as infile:
        for record in SeqIO.parse(infile, args.in_format):
            taxa_list.append(record.description)
            mat_dict[record.name] = pd.Series(list(record.seq))
    return taxa_list, mat_dict


def get_branch_lens():
    """
    Parse input tree and retrieve branch lengths
    @return:
    """
    tree = Tree(args.tree)
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
    tree = Tree(args.tree)
    nodes = [tree & x for x in mylist]
    tree.prune(nodes, preserve_branch_length=True)
    tree.write(outfile=f'{speed}.tre')

    return mylist


def trim_matrix():
    """

    @return:
    """
    with open(args.matrix, 'r') as infile, open('fast.phy', 'w') as fast_fas, open('slow.phy', 'w') as slow_fas:
        fast_recs = []
        slow_recs = []
        for record in SeqIO.parse(infile, args.in_format):
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
    i = 0
    os.mkdir(f'chunks_{args.step}')
    os.chdir(f'chunks_{args.step}')
    for chunk in range(args.step, len(sorted_sites), args.step):
        with open(f'step{i}', 'w') as res:
            records = []
            for name, seq in matrix_dict.items():
                seq = "".join(seq[sorted_sites[chunk:]].values)
                records.append(SeqRecord(Seq(seq, IUPAC.protein),
                                         id=name,
                                         name='',
                                         description=''))
            write_seqs(res, records)

        i += 1


if __name__ == '__main__':
    description = ('Removes the most heterotachous sites within a phylogenomic supermatrix in a stepwise fashion,'
                   ' leading to a user defined set of new matrices with these sites removed.')
    usage = 'heterotachy.py -t path/to/tree -m path/to/matrix [OPTIONS]'
    parser, optional, required = help_formatter.initialize_argparse(name='heterotachy.py',
                                                                    desc=description,
                                                                    usage=usage)

    required.add_argument('-t', '--tree', type=str, metavar='', required=True,
                          help=textwrap.dedent("""\
                          Path to tree
                          """))
    required.add_argument('-m', '--matrix', type=str, metavar='', required=True,
                          help=textwrap.dedent("""\
                          Path to supermatrix
                          """))
    # Optional Arguments

    optional.add_argument('-s', '--step', type=int, default=3000, metavar='N',
                          help=textwrap.dedent("""\
                          Size of removal step (i.e., 1000 sites removed) to exhaustion
                          Default: 3000
                          """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='phylip-relaxed',
                          help=textwrap.dedent("""\
                          Input format of matrix
                          Options: fasta, nexus, phylip (names truncated at 10 characters), 
                          or phylip-relaxed (names are not truncated)
                          Default: phylip-relaxed
                          """))
    optional.add_argument('-f', '--out_format', metavar='<format>', type=str, default='phylip-relaxed',
                          help=textwrap.dedent("""\
                          Desired format of the output chunks.
                          Options: fasta, nexus, phylip (names truncated at 10 characters), 
                          or phylip-relaxed (names are not truncated)
                          Default: phylip-relaxed
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

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
