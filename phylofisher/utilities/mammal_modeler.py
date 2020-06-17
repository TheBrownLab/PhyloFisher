#!/usr/bin/env python
import textwrap
import string
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from phylofisher import help_formatter


def id_generator(size=10, chars=string.ascii_uppercase):
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
    records = []
    for record in SeqIO.parse(matrix, args.in_format):
        seqs += 1
        if not length:
            length = len(record.seq)
        uname = unique_name(pseudonames)
        pseudonames[record.name] = uname
        pseudonames_rev[uname] = record.name
        records.append(SeqRecord(record.seq,
                                 id=uname,
                                 name='',
                                 description=''))
    SeqIO.write(records, 'TEMP_alph_seq.phy', 'phylip-sequential')
    return pseudonames, pseudonames_rev


def fake_tree(treefile, pseudonames):
    with open('TEMP_alph_seq.tre', 'w') as res, open('key_alph_seq.tsv', 'w') as key_tsv:
        original = open(treefile).readline()
        key_tsv.write('Original Name\tID\n')
        for key, value in pseudonames.items():
            original = original.replace(key, value)
            key_tsv.write(f'{key}\t{value}\n')
        res.write(original)


if __name__ == "__main__":
    description = 'Prepares supermatrix and tree for MAMMaL analysis'
    parser, optional, required = help_formatter.initialize_argparse(name='mammal_modeler.py',
                                                                    desc=description,
                                                                    usage='mammal_modeler.py '
                                                                          '[OPTIONS] -tr <tree_file> -m <matrix>')

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
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='phylip-relaxed',
                          help=textwrap.dedent("""\
                              Input format of matrix
                              Options: fasta, nexus, phylip (names truncated at 10 characters), 
                              or phylip-relaxed (names are not truncated)
                              Default: phylip-relaxed
                              """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

    pseudo_, pseudo_rev_ = fake_phylip(args.matrix)
    fake_tree(args.tree, pseudo_)
