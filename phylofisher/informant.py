#!/usr/bin/env python
import configparser
import glob
import os
import sys
import textwrap
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from phylofisher import help_formatter, tools


def collect_names(files):
    """
    Collect all short names. Collest information about 'route'
    for newly added sequences (BBH, SBH, HMM).
    input: fasta files with genes (from fisher)
    return: set of sorted shortnames; dictionary with genes 
    as a keys and deeper level dictionaries with shortnames of organisms 
    as keys and BBH or SBH or HMM as values
    """
    names = set()
    routes = defaultdict(dict)
    for file in files:
        gene_name = file.split('/')[-1].split('.')[0]
        for record in SeqIO.parse(file, 'fasta'):
            if record.name.count('_') >= 3:
                name = record.name.split('_')[0]
                route = record.name.split('_')[2]
                routes[gene_name][name] = route
            elif '_' in record.name:
                name = record.name.split('_')[0]
            else:
                name = record.name
            names.add(name)
    return sorted(list(names)), routes


def get_gene_column(gene, names):
    """
    Collect information about number of gene variants for all organisms.
    input: fasta file for a given gene
    return:  pd.Series with information about present(>0)/absence(0) of a gene
    for all organisms. 
    example: Albugo: 3, Naegr: 0, ...
    """
    gene_name = gene.split('/')[-1].split('.')[0]
    column = pd.Series(np.zeros(len(names)), index=names, name=gene_name)
    for record in SeqIO.parse(gene, "fasta"):
        if "_" in record.name:
            org = record.name.split('_')[0]
        else:
            org = record.name
        column[org] += 1
    return column.astype(int)


def make_table(folder):
    """
    Collect information about all genes as pd.Series and organize them
    into pd.DataFrame
    input: fodler with genes in fasta format
    return: pd.DataFrame with genes; dictionary with information about
    'routes' used for sequence selection (BBH, SBH, HMM)
    """
    genes = [gene for gene in glob.glob(f'{folder}/*.fas') if os.path.isfile(gene)]
    names, routes = collect_names(genes)
    columns = []
    for gene in genes:
        columns.append(get_gene_column(gene, names))
    df = pd.DataFrame(columns)
    df = df.transpose()
    df = df.reindex(sorted(df.columns), axis=1)  # comment me
    return df, routes


def table_with_routes(df, routes):
    """Create occupancy table and add information about routes used for sequence
    selection (BBH, SBH, HMM) to input dataframe
    input: pd.DataFrame with information about all genes (present(>0)/absent(0)) 
    for all organisms, dictionary with routes for all sequences
    output: modified input dataframe with information about route
    """
    full_names = []
    high_tax_list = []
    low_tax_list = []
    for org in in_taxa_dict.keys():
        group, subtax, long_name = in_taxa_dict[org]
        high_tax_list.append(group)
        low_tax_list.append(subtax)
        full_names.append(long_name)

    df = df[df.index.isin(in_taxa_dict.keys())]
    no_seqs = set(in_taxa_dict.keys()) - set(df.index)

    for taxon in no_seqs:
        df.loc[taxon] = len(df.columns) * [0]

    df.index.name = 'Unique ID'
    df.insert(loc=0, column='Lower Taxonomy', value=low_tax_list)
    df.insert(loc=0, column='Higher Taxonomy', value=high_tax_list)
    df.insert(loc=0, column='Full Name', value=full_names)

    df = df.sort_index(axis=0)
    df.to_csv(f'{output_fold}/occupancy.tsv', sep='\t')

    # Adds routes to df
    for gene in df.columns:
        df[gene] = df[gene].apply(str)
        for org in df[gene].index:
            if org in routes[gene]:
                df.at[org, gene] = f'{df[gene][org]}_{routes[gene][org]}'

    df.to_csv(f'{output_fold}/occupancy_with_routes.tsv', sep='\t')

    return df


def check_paralogs():
    """Collect all Unique IDs for organisms with at least one paralog
    in the dataset.
    input: None
    return: set of short names of organisms with at least one paralog
    """
    paralogs = set()
    paralog_fold = os.path.dirname(args.metadata)
    for file in glob.glob(f'{paralog_fold}/paralogs/*.fas'):
        for record in SeqIO.parse(file, 'fasta'):
            paralogs.add(record.name.split('.')[0])
    return paralogs


def get_routes():
    """

    :return:
    """
    my_routes = dict()
    for org in in_taxa_dict.keys():
        sbh = {'SBH': 0, 'BBH': 0, 'HMM': 0}
        for val in res.loc[org].values[2:]:
            if '_' in val:
                count, route = val.split('_')
                sbh[route] += int(count)

        my_routes[org] = [sbh['SBH'], sbh['BBH'], sbh['HMM']]

    return my_routes


def stats_orgs(df, new_data=False):
    """
    Create tsv file with basic summary about analyzed dataset without information
    about 'routes' and paralogs.
    input: dataframe
    return: None
    """
    rows = []

    if new_data:
        df = df[df.index.isin(in_taxa_dict.keys())]
    else:
        df = df[df.index.isin(db_taxa_dict.keys())]

    df2 = df.copy()
    df2[df2 >= 1] = 1

    df = df.sum(axis=1).to_frame()

    if new_data:
        df[f"Genes out of {len(matrix.columns)}"] = df2.sum(axis=1).to_frame()
        df = df.rename(columns={0: f"Sequences Collected"})

    else:
        df = df.rename(columns={0: f"Genes out of {len(matrix.columns)}"})

    # Fill in taxonomic information
    if new_data:
        list_of_dicts = [{key: value[i] for key, value in in_taxa_dict.items()} for i in range(3)]
    else:
        list_of_dicts = [{key: value[i] for key, value in db_taxa_dict.items()} for i in range(3)]
    df['Long Name'] = df.index.map(list_of_dicts[2])
    df['Higher Taxonomy'] = df.index.map(list_of_dicts[0])
    df['Lower Taxonomy'] = df.index.map(list_of_dicts[1])

    # Rearrange Columns to Put Genes after taxa stats
    cols = df.columns.tolist()
    cols = cols[2:] + cols[:2]
    df = df[cols]

    if new_data:
        routes_dict = get_routes()
        list_of_routes_dicts = [{key: value[i] for key, value in routes_dict.items()} for i in range(3)]
        df["#SBH"] = df.index.map(list_of_routes_dicts[0])
        df["#BBH"] = df.index.map(list_of_routes_dicts[1])
        df["#HMM"] = df.index.map(list_of_routes_dicts[2])
        out_filename = 'new_taxa_stats.tsv'
    else:
        out_filename = 'db_taxa_stats.tsv'

    # Fill in columns for including in SGT construction. By default all are yes
    has_paralogs = check_paralogs()
    if new_data:
        sgt_dict = {org: 'yes' for org in in_taxa_dict.keys()}
    else:
        sgt_dict = {org: 'yes' for org in db_taxa_dict.keys()}
    df['SGT'] = df.index.map(sgt_dict)

    # Fill in column for paralogs. If no paralogs entry is 'none'.
    # If there are paralogs entry is 'yes'. If there are paralogs, but --ortholog_only is given entry is 'no'.
    if new_data:
        pass
    else:
        paralogs_dict = {org: ('yes' if org in has_paralogs and not args.orthologs_only
                               else 'no' if org in has_paralogs and args.orthologs_only else 'none')
                         for org in db_taxa_dict}
        df['Paralogs'] = df.index.map(paralogs_dict)

    df = df.rename_axis('Unique ID')
    df.to_csv(f'{output_fold}/{out_filename}', sep='\t')


def stats_gene(df):
    """Write basic summary for all genes to a file
    input: dataframe
    return: None
    """
    # res.write(f'Gene,Total,total[%],SGT\n')
    taxa_count = len(df)
    df = df.sum().to_frame()
    df = df.rename(columns={0: 'Number of Taxa'})
    df[f'Percent of Total Taxa (out of {taxa_count})'] = round((df['Number of Taxa'] / taxa_count) * 100, 2)
    df = df.rename_axis('Gene Name')
    df = df.sort_values(by=['Number of Taxa'], ascending=False)
    df['SGT'] = ['yes'] * len(df)
    df.to_csv(f'{output_fold}/gene_stats.tsv', sep='\t')


if __name__ == '__main__':
    description = 'Produce preliminary statistics about newly input data.'
    parser, optional, required = help_formatter.initialize_argparse(name='informant.py',
                                                                    desc=description,
                                                                    usage='informant.py -i input_folder [OPTIONS]')

    # Optional Arguments
    optional.add_argument('--orthologs_only', action='store_true',
                          help=textwrap.dedent("""\
                          Paralogs will NOT be included from any taxa in the starting 
                          database in downstream single gene tree construction.
                          """))

    in_help = 'Path to fisher.py output directory'
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, in_help=in_help, out_dir=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    args.input_metadata = os.path.abspath(config['PATHS']['input_file'])
    args.metadata = os.path.join(dfo, 'metadata.tsv')

    output_fold = os.path.basename(os.path.normpath(args.input)) + '/informant_stats'
    if os.path.isdir(output_fold):
        sys.exit(f'Error: {output_fold} folder already exists.')
    else:
        os.mkdir(output_fold)

    db_taxa_dict = tools.parse_metadata(args.metadata)
    in_taxa_dict = tools.parse_metadata(args.input_metadata, input_meta=True)
    all_taxa_dict = {**db_taxa_dict, **in_taxa_dict}

    matrix, routes = make_table(args.input)

    res = table_with_routes(matrix, routes)

    stats_orgs(matrix)
    stats_orgs(matrix, new_data=True)
    stats_gene(matrix)
