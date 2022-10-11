#!/usr/bin/env python

'''
Calculates the observed occurrences of clades of interest in bootstrap trees.
'''

import configparser
import os
import shutil
import sys
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from ete3 import Tree

from phylofisher import help_formatter

plt.rcParams["figure.figsize"] = (10, 5)


def bipartitions(tree):
    '''
    Gets the bipartitions of a tree

    :param tree: tree
    :type tree: ETE3 Tree object
    :return: bipartitions
    :rtype: set
    '''
    bipar = set()
    all_ = frozenset(tree.get_leaf_names())
    for node in tree.traverse('preorder'):
        if node.is_leaf() is False:
            taxa = frozenset(node.get_leaf_names())
            rest = frozenset(taxa.symmetric_difference(all_))
            if len(taxa) > 1:
                bipar.add(taxa)
            if len(rest) > 1:
                bipar.add(rest)
    return bipar


def support(trees):
    '''
    Get support from trees

    :param trees: path to trees
    :type trees: str
    :return: support dictionary, all taxa
    :rtype: dict, list
    '''
    bootstrap = []
    n_trees = 0
    all_taxa = set()
    with open(trees, 'r', encoding='utf8') as tree_file:
        for tree in tree_file:
            tree = tree.strip()
            leaves = list(Tree(tree).get_leaf_names())
            all_taxa.update(leaves)

    with open(trees, 'r', encoding='utf8') as tree_file:
        for tree in tree_file:
            tree = tree.strip()
            tree = Tree(tree)
            n_trees += 1
            bootstrap = bootstrap + list(bipartitions(tree))
    norm_bootstrap = {}
    counted = dict(Counter(bootstrap))
    for key, value in counted.items():
        norm_bootstrap[key] = value / n_trees

    return norm_bootstrap, all_taxa


def get_support(group, supp_dict):
    '''
    Checks support for a group

    :param group: group of interest
    :type group: set
    :param supp_dict: support dictionary
    :type supp_dict: dict
    :return: support
    :rtype: float
    '''
    query = frozenset(group)
    if query in supp_dict:
        return supp_dict[query]

    return 0


def get_taxa_in_group(groups):
    '''
    Gets taxa in a group

    :param groups: group of interest
    :type groups: list
    :return: taxa in group
    :rtype: list
    '''
    df = pd.read_csv(METADATA, delimiter='\t')
    taxa = []
    for group in groups:
        if group in list(df['Higher Taxonomy']):
            taxa = taxa + list(df.loc[df['Higher Taxonomy'] == group, 'Unique ID'])
        elif group in list(df['Lower Taxonomy']):
            taxa = taxa + list(df.loc[df['Lower Taxonomy'] == group, 'Unique ID'])
        else:
            sys.exit(f'{group} is not in the database\'s metadata')

        if args.chimeras:
            for chim_key, chim_value in chim_dict.items():
                if group in (chim_value[0], chim_value[1]):
                    taxa.append(chim_key)

    return taxa


def parse_groups(input_file):
    '''
    Parse groups file

    :param input_file: path to groups file
    :type input_file: str
    :return: groups
    :rtype: dict
    '''
    query_dict = {}
    with open(input_file, 'r', encoding='utf8') as group_file:
        for line in group_file:
            line = line.strip()
            if ':' in line:
                group = line.split(':')[0]
                orgs = line.split(':')[1]
                query_dict[group] = [org.strip() for org in orgs.split(',')]
            else:
                groups = list(line.split('+'))
                query_dict[line] = get_taxa_in_group(groups)
    return query_dict.items()


def file_to_series(file):
    '''
    Converts a file to a pandas series

    :param file: path to input file
    :type file: str
    :return: pandas series
    :rtype: pd.Series
    '''
    print(file)
    sup_dict, all_taxa = support(file)
    group_sup = {}
    for group, orgs in queries:
        orgs = set(orgs) - ((set(orgs) - all_taxa))
        group_sup[group] = get_support(orgs, sup_dict)
    s = pd.Series(group_sup)
    return s


def parse_bss():
    '''
    Parse bss files

    :return: rows of bss file
    :rtype: list
    '''
    rows = []
    count = 0
    with open(args.bs_files, 'r', encoding='utf8') as bs_file:
        for row in bs_file:
            row = row.strip()
            if not os.path.isfile(row):
                raise FileNotFoundError(f'{row} does not exist')
            row_series = file_to_series(row.strip())
            row_series.name = os.path.basename(row)
            # row.name = f'Step {n}'
            rows.append(row_series)
            count += 1

    return rows


def main():
    '''
    Makes plots of bipartition support

    :return: bipartition dataframe
    :rtype: pandas dataframe
    '''
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
    fig.savefig(f'{args.output}/bipartition_examiner.pdf')
    df.to_csv(f'{args.output}/bipartition_examiner.tsv', sep='\t')
    return df


if __name__ == "__main__":
    DESCRIPTION = 'Calculates the observed occurrences of clades of interest in bootstrap trees.'
    parser, optional, required = help_formatter.initialize_argparse(
        name='bipartition_examiner.py',
        desc=DESCRIPTION,
        usage='bipartition_examiner.py [OPTIONS] -i /path/to/input/'
        )

    # Required Arguments
    required.add_argument('-b', '--bs_files', required=True, type=str, metavar='<bs_files>',
                          help=textwrap.dedent("""\
                          Path bootstrap files.
                          """))
    required.add_argument('-g', '--groups', type=str, required=True, metavar='<groups.txt>',
                          help=textwrap.dedent("""\
                          Path to text file containing groups of interest.
                          """))

    # Optional Arguments
    optional.add_argument('--database', type=str, metavar='<path/to/db>',
                          help=textwrap.dedent("""\
                          Path to database if not using config.ini
                          """))
    optional.add_argument('--no_db', action='store_true',
                          help=textwrap.dedent("""\
                          Do not use a data base
                          """))
    optional.add_argument('--chimeras', type=str, metavar='<path/to/chimeras>',
                          help=textwrap.dedent("""\
                          A .tsv containing a Unique ID, higher taxonomy, and lower taxonomy for each chimera within the input bootstrap files.

                          Example:
                            Chimera_ID 1\tHigher_Tax\tLower_Tax
                            Chimera_ID 2\tHigher_Tax\tLower_Tax
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
        elif os.path.isfile('config.ini'):
            config = configparser.ConfigParser()
            config.read('config.ini')
            dfo = str(Path(config['PATHS']['database_folder']).resolve())
        METADATA = str(os.path.join(dfo, 'metadata.tsv'))

    # Make output directory
    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    # Parse chimera file
    if args.chimeras:
        chim_dict = {}
        with open(args.chimeras, 'r', encoding='utf8') as infile:
            for line in infile:
                line = line.strip()
                split_line = line.split('\t')
                chim_dict[split_line[0]] = split_line[1:3]

    queries = parse_groups(args.groups)
    main()
