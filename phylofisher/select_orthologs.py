#!/usr/bin/env python
import configparser
import os
import textwrap
from pathlib import Path

from phylofisher import help_formatter, tools


def in_out_completeness(df):
    """

    :param df:
    :return:
    """
    # Gets all taxa unique IDS
    all_taxa_set = set()
    with open(f'{dfo}/metadata.tsv', 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')[0]
            all_taxa_set.add(line)

    # Gets out group Unique IDS
    out_group_set = set()
    with open(args.out_group, 'r') as infile:
        for line in infile:
            line = line.strip()
            out_group_set.add(line)

    # Determine in group from all taxa - out group
    in_group_set = all_taxa_set - out_group_set

    # Calculates completeness for each gene based on in group presence
    in_df = df[df.index.isin(list(in_group_set))]
    in_df = in_df.sum().divide(other=len(in_df))

    # Calculates completeness for each gene based on out group presence
    out_df = df[df.index.isin(list(out_group_set))]
    out_df = out_df.sum().divide(other=len(out_df))

    return out_df.to_dict(), in_df.to_dict()


def make_subset_tsv():
    """

    :return: df
    """
    df = gene_comp.to_frame()
    df = df.rename(columns={0: 'Completeness'})
    df = df.sort_values(by=['Completeness'], ascending=False)
    df = df.round({'Completeness': 3})

    if args.out_group:
        out_dict, in_dict = in_out_completeness(matrix)
        df['In-Group Completeness'] = df.index.map(in_dict)
        df['Out-Group Completeness'] = df.index.map(out_dict)
        df = df.round({'In-Group Completeness': 3})
        df = df.round({'Out-Group Completeness': 3})

    to_keep = ['yes' for k in df.index]
    df['Include in Subset'] = to_keep

    if args.gene_number:
        num_to_keep = args.gene_number
        i = 1
        for gene, row in df.iterrows():
            if i > num_to_keep:
                df.at[gene, 'Include in Subset'] = 'no'
            i += 1

    elif args.percent_complete:
        if args.percent_complete > 1:
            args.percent_complete = args.percent_complete / 100

        if args.out_group:
            for gene, row in df.iterrows():
                if df.at[gene, 'In-Group Completeness'] < args.percent_complete:
                    df.at[gene, 'Include in Subset'] = 'no'
        else:
            for gene, row in df.iterrows():
                if df.at[gene, 'Completeness'] < args.percent_complete:
                    df.at[gene, 'Include in Subset'] = 'no'

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
    optional.add_argument('--out_group', type=str, metavar='out_group.txt', default=None,
                          help=textwrap.dedent("""\
                          Path to text file containing out-group taxa Unique IDs.
                          """))
    optional.add_argument('-n', '--gene_number', type=int, metavar='<N>', default=None,
                          help=textwrap.dedent("""\
                          Number of genes in subset.
                          Cannot be used with percent_complete.
                          """))
    optional.add_argument('-c', '--percent_complete', type=float, metavar='<N>', default=None,
                          help=textwrap.dedent("""\
                          Threshold for percent complete when subsetting.
                          Cannot be used with gene_number.
                          """))
    optional.add_argument('--chimeras', type=str, metavar='chimeras.tsv', default=None,
                          help=textwrap.dedent("""\
                          A TSV containing chimeras, higher and lower taxonomic designations, 
                          and the taxa comprising each chimera.
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    metadata = f'{dfo}/metadata.tsv'
    orthologs_dir = f'{dfo}/orthologs/'

    matrix = tools.completeness(args=args, input_dir=orthologs_dir, genes=True)
    if os.path.isfile('select_taxa.tsv'):
        matrix = update_df_taxa(matrix)
    taxa_count, _ = matrix.shape
    gene_comp = matrix.sum().divide(other=taxa_count)

    make_subset_tsv()

    tools.make_plot(gene_comp, f'gene_comp', y_count=taxa_count, genes=True)
