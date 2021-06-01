#!/usr/bin/env python
import configparser
import os
import textwrap
import pandas as pd
from pathlib import Path

from phylofisher import help_formatter, tools


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
    taxa = tools.parse_metadata(metadata)
    df = taxa_comp.to_frame()
    df = df.rename(columns={0: 'Completeness'})

    if args.chimeras:
        with open(args.chimeras, 'r') as infile:
            for line in infile:
                line = line.strip()
                split_line = line.split('\t')
                taxa[split_line[0]] = split_line[1:3]

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


def update_df_ortho(df):
    """

    :return:
    """
    with open('select_orthologs.tsv', 'r') as infile:
        infile.readline()
        to_drop = []
        for line in infile:
            line = line.strip()
            gene, _, include = line.split()
            if include == 'no':
                to_drop.append(gene)

    df = df.drop(to_drop)
    return df


def gen_chimera(df):
    with open(args.chimeras, 'r') as infile:
        chim_dict = {}
        for line in infile:
            line = line.strip()
            split_line = line.split('\t')
            chim_dict[split_line[0]] = split_line[3:]

    for key in chim_dict.keys():
        df['int_chim'] = 0
        for taxon in chim_dict[key]:
            df['int_chim'] = df['int_chim'] + df[taxon]
            df = df.drop(taxon, axis=1)

        df.loc[df['int_chim'] >= 1, key] = 1
        df.loc[df['int_chim'] < 1, key] = 0
        df = df.drop('int_chim', axis=1)

    return df


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
    optional.add_argument('--to_exclude', type=str, metavar='exc_taxa.txt', default=None,
                          help=textwrap.dedent("""\
                          List of taxa to exclude.
                          """))
    optional.add_argument('--to_include', type=str, metavar='inc_taxa.txt', default=None,
                          help=textwrap.dedent("""\
                          List of taxa to include.
                          """))
    optional.add_argument('--chimeras', type=str, metavar='chimeras.tsv', default=None,
                          help=textwrap.dedent("""\
                           A .tsv containing a Unique ID, higher and lower taxonomic designations, 
                           and the Unique IDs of the taxa to collapse, for each chimera one per line
                          """))

    args = help_formatter.get_args(parser, optional, required, inp_dir=False, pre_suf=False, out_dir=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    metadata = f'{dfo}/metadata.tsv'
    orthologs_dir = f'{dfo}/orthologs/'

    matrix = tools.completeness(args=args, input_dir=orthologs_dir, genes=True)
    matrix = matrix.transpose()

    if os.path.isfile('select_orthologs.tsv'):
        matrix = update_df_ortho(matrix)

    if args.chimeras:
        matrix = gen_chimera(matrix)
        matrix.to_csv('text_df.csv')

    gene_count, _ = matrix.shape
    taxa_comp = matrix.sum().divide(other=gene_count)

    make_subset_tsv()
    tools.make_plot(taxa_comp, f'taxa_comp', y_count=gene_count, genes=False)
