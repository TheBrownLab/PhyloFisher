#!/usr/bin/env python

import configparser
import os
import textwrap
from pathlib import Path

from phylofisher import help_formatter, subset_tools


def in_out_completeness(df):
    """

    :param df:
    :return:
    """
    pass


def make_subset_tsv():
    df = gene_comp.to_frame()
    df = df.rename(columns={0: 'Completeness'})
    df = df.sort_values(by=['Completeness'], ascending=False)
    df = df.round({'Completeness': 3})
    to_keep = ['no' for k in df.index]
    df['Include in Subset'] = to_keep

    if args.gene_number:
        num_to_keep = args.gene_number

    elif args.percent_complete:
        if args.percent_complete > 1:
            args.percent_complete = args.percent_complete / 100
        num_to_keep = int(round(args.percent_complete * len(df)))

    else:
        num_to_keep = len(df)

    i = 1
    for gene, row in df.iterrows():
        if i <= num_to_keep:
            df.at[gene, 'Include in Subset'] = 'yes'
        i += 1

    df.to_csv(f'select_orthologs.tsv', sep='\t')
    return df


def update_df_taxa(df):
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
    description = 'Selects orthologs to be included in super matrix construction'
    parser, optional, required = help_formatter.initialize_argparse(name='select_orthologs.py',
                                                                    desc=description,
                                                                    usage='select_orthologs.py '
                                                                          '[OPTIONS]')

    # Optional Arguments
    optional.add_argument('-n', '--gene_number', type=int, metavar='<N>', default=None,
                          help=textwrap.dedent("""\
                          Number of genes in subset.
                          This will be ignored if not used with --subset.
                          Cannot be used with percent_complete.
                          """))
    optional.add_argument('-c', '--percent_complete', type=float, metavar='<N>', default=None,
                          help=textwrap.dedent("""\
                          Threshold for percent missing when subsetting.
                          This will be ignored if not used with --subset.
                          Cannot be used with gene_number.
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    metadata = f'{dfo}/metadata.tsv'
    orthologs_dir = f'{dfo}/orthologs/'

    matrix = subset_tools.completeness(args=args, input_dir=orthologs_dir, genes=True)
    if os.path.isfile('select_taxa.tsv'):
        matrix = update_df_taxa(matrix)
    taxa_count, _ = matrix.shape
    gene_comp = matrix.sum().divide(other=taxa_count)

    make_subset_tsv()

    subset_tools.make_plot(gene_comp, f'gene_comp', y_count=taxa_count, genes=True)
