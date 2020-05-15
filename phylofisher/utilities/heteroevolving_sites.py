#!/usr/bin/env python
import os
import textwrap
from ete3 import Tree
from Bio import SeqIO

from phylofisher import help_formatter
from phylofisher.utilities.fast_tax_removal import Leaves


def fast_and_slow(matrix, portion, sorted_orgs):
    amount = round(len(sorted_orgs)*portion)
    slow = set(sorted_orgs[(len(sorted_orgs)-amount):])
    fast = set(sorted_orgs[:amount])
    print(len(sorted_orgs))
    print(len(slow))
    print(len(fast))
    with open(f'{args.output}/slow.fas', 'w') as s, open(f'{args.output}/fast.fas', 'w') as f:
        slow_seqs, fast_seq = [], []
        for record in SeqIO.parse(matrix, args.in_format):
            if record.name in slow:
                s.write(f'>{record.name}\n{record.seq}\n')
            elif record.name in fast:
                f.write(f'>{record.name}\n{record.seq}\n')


if __name__ == '__main__':
    description = 'Fast taxon removal'
    parser, optional, required = help_formatter.initialize_argparse(name='heteroevolving_sites.py',
                                                                    desc=description,
                                                                    usage='heteroevolving_sites.py [OPTIONS] '
                                                                          '-t <tree> '
                                                                          '-m <matrix>')
    # Required Arguments
    required.add_argument('-m', '--matrix', required=True, type=str, metavar='matrix',
                          help=textwrap.dedent("""\
                              Path to input matrix.
                              """))
    required.add_argument('-t', '--tree', required=True, type=str, metavar='tree',
                          help=textwrap.dedent("""\
                              Path to input tree.
                              """))

    # Optional Arguments
    optional.add_argument('-in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                              Input matrix format if not FASTA.
                              Options: fasta, phylip (names truncated at 10 characters), 
                              phylip-relaxed (names are not truncated), or nexus.
                              Default: fasta
                              """))
    # optional.add_argument('-out_format', metavar='<format>', type=str, default='fasta',
    #                       help=textwrap.dedent("""\
    #                                   Desired output format.
    #                                   Options: fasta, phylip (names truncated at 10 characters),
    #                                   phylip-relaxed (names are not truncated), or nexus.
    #                                   Default: fasta
    #                                   """))

    optional.add_argument('-p', '--portion', default=0.2, type=float,
                          help=textwrap.dedent("""\
                              Portion
                              Default: 0.2
                              """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

    os.mkdir(args.output)
    tree = Tree(args.tree)
    leaves = Leaves(tree)
    sorted_taxa = leaves.org_speed()
    fast_and_slow(args.matrix, args.portion, sorted_taxa)