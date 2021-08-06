#!/usr/bin/env python
import configparser
import os
import sys
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from ete3 import Tree

from phylofisher import help_formatter

plt.rcParams["figure.figsize"] = (10, 5)


def get_taxa_set(trees):
    """
    """


def bipartitions(tree):
    """

    :param tree:
    :return:
    """
    bipar = set()
    all_ = frozenset(tree.get_leaf_names())
    for node in tree.traverse('preorder'):
        if node.is_leaf() is False:
            taxa = frozenset(node.get_leaf_names())
            if len(taxa) < (len(all_) - len(taxa)):
                if len(taxa) > 1:
                    bipar.add(taxa)
            else:
                rest = taxa.symmetric_difference(all_)
                if len(rest) > 1:
                    bipar.add(rest)
    return bipar


def support(trees):
    """

    :param trees:
    :return:
    """
    bootstrap = []
    n_trees = 0
    all_taxa = set()
    for line in open(trees):
        line = line.strip()
        leaves = list(Tree(line).get_leaf_names())
        all_taxa.update(leaves)

    for line in open(trees):
        line = line.strip('')
        tree = Tree(line)
        n_trees += 1
        bootstrap = bootstrap + list(bipartitions(tree))
    norm_bootstrap = {}
    counted = dict(Counter(bootstrap))
    for key, value in counted.items():
        norm_bootstrap[key] = value / n_trees
    
    return norm_bootstrap, all_taxa


def get_support(group, supp_dict):
    """

    :param group:
    :param supp_dict:
    :return:
    """
    query = frozenset(group)
    if query in supp_dict:
        return supp_dict[query]
    else:
        return 0


def get_taxa_in_group(groups):
    """

    :return:
    """
    df = pd.read_csv(metadata, delimiter='\t')
    taxa = []
    for group in groups:
        if group in list(df['Higher Taxonomy']):
            taxa = taxa + list(df.loc[df['Higher Taxonomy'] == group, 'Unique ID'])
        elif group in list(df['Lower Taxonomy']):
            taxa = taxa + list(df.loc[df['Lower Taxonomy'] == group, 'Unique ID'])
        else:
            sys.exit(f'{group} is not in the database\'s metadata')

        if args.chimeras:
            for chimera in chim_dict.keys():
                if chim_dict[chimera][0] == group or chim_dict[chimera][1] == group:
                    taxa.append(chimera)

    return taxa


def parse_groups(input_file):
    """

    :param input_file:
    :return:
    """
    query_dict = {}
    for line in open(input_file):
        line = line.strip()
        if ':' in line:
            group = line.split(':')[0]
            orgs = line.split(':')[1]
            query_dict[group] = [org.strip() for org in orgs.split(',')]
        else:
            group = line
            groups = [group for group in group.split('+')]
            query_dict[group] = get_taxa_in_group(groups)
    return query_dict.items()


def file_to_series(file):
    """

    :param file:
    :return:
    """
    print(file)
    sup_dict, all_taxa = support(file)
    group_sup = {}
    for group, orgs in queries:
        orgs = set(orgs) - ((set(orgs) - all_taxa))
        group_sup[group] = get_support(orgs, sup_dict)
    s = pd.Series(group_sup)
    return s


def parse_bss():
    """

    :return:
    """
    columns = []
    n = 0
    for line in open(args.bs_files):
        line = line.strip()
        column = file_to_series(line.strip())
        column.name = os.path.basename(line)
        # column.name = f'Step {n}'
        columns.append(column)
        n += 1
        
    return columns


def main():
    """

    :return:
    """
    columns = parse_bss()
    df = pd.DataFrame(columns)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if args.bar_plot:
        df = df.transpose()
        df.plot.bar(figsize=(20, 10),
                    legend=True,
                    fontsize=12,
                    ax=ax,
                    cmap="Set2")
        ax.set_xlabel('Groups', size=16)

    else:
        df.plot(figsize=(35, 10),
                legend=True,
                fontsize=8,
                ax=ax,
                cmap="Set2")
        ax.set_xticks(range(len(list(df.index))))

    ax.set_xticklabels(list(df.index), rotation=45, ha='right')
    ax.set_ylabel('Bootstrap Support', size=12)
    ax.legend(loc=0, prop={'size': 12})
    ax.margins(x=0.005)
    plt.tight_layout()
    fig.savefig("bipartition_examiner.pdf")
    df.to_csv("bipartition_examiner.tsv", sep='\t')
    return df


if __name__ == "__main__":
    description = 'Calculates the observed occurrences of clades of interest in bootstrap trees.'
    parser, optional, required = help_formatter.initialize_argparse(name='bipartition_examiner.py',
                                                                    desc=description,
                                                                    usage='bipartition_examiner.py '
                                                                          '[OPTIONS] -i /path/to/input/')

    # Required Arguments
    required.add_argument('-b', '--bs_files', required=True, type=str, metavar='',
                          help=textwrap.dedent("""\
                          Path bootstrap files.
                          """))
    required.add_argument('-g', '--groups', type=str, required=True, metavar='',
                          help=textwrap.dedent("""\
                          groups
                          """))

    # Optional Arguments
    optional.add_argument('-on', '--output_name', type=str, metavar='out', default='out',
                          help=textwrap.dedent("""\
                          """))
    optional.add_argument('--database', type=str, metavar='path/to/db',
                          help=textwrap.dedent("""\
                          Path to database if not using config.ini
                          """))
    optional.add_argument('--no_db', action='store_true',
                          help=textwrap.dedent("""\
                          
                          """))
    optional.add_argument('--chimeras', type=str, metavar='chimera_info.tsv',
                          help=textwrap.dedent("""\
                          Path to TSV file containing taxonomic designations for chimeras.
                          """))
    optional.add_argument('--bar_plot', action='store_true',
                          help=textwrap.dedent("""\
                          Plot categorical data as a barplot.
                          Default: Plot series data as a line graph.
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

    if not args.no_db:
        if args.database:
            dfo = os.path.abspath(args.database)
        else:
            config = configparser.ConfigParser()
            config.read('config.ini')
            dfo = str(Path(config['PATHS']['database_folder']).resolve())

        metadata = str(os.path.join(dfo, 'metadata.tsv'))

    if args.chimeras:
        chim_dict = {}
        with open(args.chimeras, 'r') as infile:
            for line in infile:
                line = line.strip()
                split_line = line.split('\t')
                chim_dict[split_line[0]] = split_line[1:3]

    queries = parse_groups(args.groups)
    main()
