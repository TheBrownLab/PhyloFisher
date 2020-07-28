#!/usr/bin/env python

import configparser
import os
import textwrap
from pathlib import Path
import pandas as pd

from phylofisher import help_formatter, subset_tools


def parse_user_inc_exc(input_file):
    """

    :return:
    """
    with open(input_file, 'r') as infile:
        user_set = set()
        for line in infile:
            line = line.strip()
            user_set.add(line)

        return user_set


def make_subset_tsv():
    """

    :return:
    """
    taxa = subset_tools.parse_metadata(metadata)
    df = taxa_comp.to_frame()
    df = df.rename(columns={0: 'Completeness'})

    # Add Taxonomic Groups to DataFrame
    high_tax_list = [taxa[ind][0] for ind in df.index]
    df['Higher Taxonomy'] = high_tax_list
    low_tax_list = [taxa[ind][1] for ind in df.index]
    df['Lower Taxonomy'] = low_tax_list

    # Reorder DataFrame putting completeness in last column and sort
    cols = df.columns.tolist()
    cols = cols[1:] + cols[:1]
    df['Completeness'] = df['Completeness'] * 100
    df['Completeness'] = df['Completeness'].round(2)
    df = df[cols]
    df = df.sort_values(by=['Higher Taxonomy', 'Lower Taxonomy', 'Completeness'])

    # Add "Include in Subset" column to DataFrame and pre-fill out
    to_keep = ['yes' for k in df.index]
    df['Include in Subset'] = to_keep

    if args.to_include:
        include_set = parse_user_inc_exc(args.to_include)
    else:
        include_set = set()

    if args.to_exclude:
        exclude_set = parse_user_inc_exc(args.to_exclude)
    else:
        exclude_set = set()

    for taxon, row in df.iterrows():
        if (row['Lower Taxonomy'] in exclude_set) or (row['Higher Taxonomy'] in exclude_set) or (taxon in exclude_set):
            df.at[taxon, 'Include in Subset'] = 'no'

        if (row['Lower Taxonomy'] in include_set) or (row['Higher Taxonomy'] in include_set) or (taxon in include_set):
            df.at[taxon, 'Include in Subset'] = 'yes'

    df.to_csv(f'select_taxa.tsv', sep='\t')


def update_dataframe(df):
    """

    :return:
    """
    with open('select_taxa.tsv', 'r') as infile:
        infile.readline()
        to_drop = []
        for line in infile:
            line = line.strip()
            taxon, _, _, _, include = line.split()
            if include == 'no':
                to_drop.append(taxon)

    df = df.drop(to_drop)

    return df


if __name__ == '__main__':
    description = 'Selects taxa to be included in super matrix construction'
    parser, optional, required = help_formatter.initialize_argparse(name='select_taxa.py',
                                                                    desc=description,
                                                                    usage='select_taxa.py '
                                                                          '[OPTIONS]')

    # Optional Arguments
    optional.add_argument('--to_include', type=str, metavar='inc_taxa.txt', default=None,
                          help=textwrap.dedent("""\
                          List of taxa to include.
                          """))
    optional.add_argument('--to_exclude', type=str, metavar='exc_taxa.txt', default=None,
                          help=textwrap.dedent("""\
                          List of taxa to exclude
                          """))

    args = help_formatter.get_args(parser, optional, required, inp_dir=False, pre_suf=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    metadata = f'{dfo}/metadata.tsv'
    orthologs_dir = f'{dfo}/orthologs/'

    taxa_comp, gene_count = subset_tools.completeness(args=args, input_dir=orthologs_dir, genes=False)
    make_subset_tsv()
    taxa_comp = update_dataframe(taxa_comp)
    subset_tools.make_plot(taxa_comp, f'taxa_comp', y_count=gene_count, genes=False)
