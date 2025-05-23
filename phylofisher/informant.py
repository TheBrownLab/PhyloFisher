#!/usr/bin/env python
import configparser
import glob
import os
import sys
import textwrap
from collections import defaultdict
from pathlib import Path
from peewee import *
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from Bio import SeqIO
from phylofisher.db_map import database, BaseModel, Genes, Taxonomies, Metadata, Sequences


from phylofisher import help_formatter, tools


def collect_names(files):
    '''
    Collect all short names. Collest information about 'route' for newly added sequences (BBH, SBH, HMM).

    :param files: list of fasta files
    :type files: list
    :return: set of sorted shortnames; dictionary with genes
    :rtype: tuple(set, dict)
    '''
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
    '''
    Collect information about number of gene variants for all organisms.

    :param gene: path to fasta file for a given gene
    :type gene: str
    :param names: set of short names of organisms
    :type names: set
    :return: pd.Series with information about present(>0)/absence(0) of a gene
    :rtype: pd.Series
    '''
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
    '''
    Collect information about all genes as pd.Series and organize them into pd.DataFrame

    :param folder: directory with fasta files for all genes
    :type folder: str
    :return: pd.DataFrame with genes; dictionary with information about 'routes' used for sequence selection (BBH, SBH, HMM)
    :rtype: tuple(pd.DataFrame, dict)
    '''
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
    '''
    Create occupancy table and add information about routes used for sequence selection (BBH, SBH, HMM) to input dataframe

    :param df: pd.DataFrame with information about all genes (present(>0)/absent(0))
    :type df: pd.DataFrame
    :param routes: dictionary with information about 'routes' used for sequence selection (BBH, SBH, HMM)
    :type routes: dict
    :return: pd.DataFrame with genes and information about routes
    :rtype: pd.DataFrame
    '''
    # Convert to boolean matrix
    df = df[df.index.isin(in_taxa_dict.keys())]

    # Get taxa with no sequences and fill in with zerox
    no_seqs = set(in_taxa_dict.keys()) - set(df.index)
    for taxon in no_seqs:
        df.loc[taxon] = len(df.columns) * [0]

    # Set index to Unique ID
    df.index.name = 'Unique ID'
    
    # Initialize new columns
    df.insert(loc=0, column='Lower Taxonomy', value=['NA'] * len(df))
    df.insert(loc=0, column='Higher Taxonomy', value=['NA'] * len(df))
    df.insert(loc=0, column='Full Name', value=['NA'] * len(df))

    for k, v in in_taxa_dict.items():
        df.at[k, 'Full Name'] = v[2]
        df.at[k, 'Higher Taxonomy'] = v[1]
        df.at[k, 'Lower Taxonomy'] = v[0]

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
    '''
    Collect all Unique IDs for organisms with at least one paralog in the dataset.

    :return: set of short names of organisms with at least one paralog
    :rtype: set
    '''
    paralogs = set()
    database.init(os.path.join(dfo, 'phylofisher.db'))
    database.connect()
    db_query = Sequences.select(Sequences.metadata).where(Sequences.is_paralog == True)
    for result in db_query:
        org = Metadata.get(Metadata.id == result.metadata).short_name
        paralogs.add(org)
    database.close()
    return paralogs


def get_routes():
    '''
    Collect information about 'routes' used for sequence selection (BBH, SBH, HMM) for all organisms.

    :return: dictionary with information about 'routes' used for sequence selection (BBH, SBH, HMM)
    :rtype: dict
    '''
    my_routes = dict()
    for org in in_taxa_dict.keys():
        sbh = {'SBH': 0, 'BBH': 0, 'HMM': 0}
        for val in res.loc[org].values[2:]:
            if '_' in val:
                count, route = val.split('_')
                sbh[route] += int(count)

        my_routes[org] = [sbh['SBH'], sbh['BBH'], sbh['HMM']]

    return my_routes


def update_homolog_tree(df):
    '''
    Update homolog tree column to only include those organisms that are provided by the user.

    :param df: pd.DataFrame with information about all genes and organisms
    :type df: pd.DataFrame
    :return: pd.DataFrame with updated homolog tree information
    :rtype: pd.DataFrame
    '''
    ht_include = set()
    with open(args.ht_include, 'r') as f:
        for line in f:
            ht_include.add(line.strip())
    
    for index, _ in df.iterrows():
        if index not in ht_include:
            df.at[index, 'Homolog Tree'] = 'no'
            df.at[index, 'Paralogs'] = 'no'
         
    return df


def stats_orgs(df, new_data=False):
    '''
    Create tsv file with basic summary about analyzed dataset without information about 'routes' and paralogs.

    :param df: pd.DataFrame with information about all genes and organisms
    :type df: pd.DataFrame
    :param new_data: Boolean indicating if the data is new or from the database
    :type new_data: bool, optional
    '''
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

    # Fill in columns for including in homolog tree construction. By default all are yes
    has_paralogs = check_paralogs()
    if new_data:
        homolog_tree_dict = {org: 'yes' for org in in_taxa_dict.keys()}
    else:
        homolog_tree_dict = {org: 'yes' for org in db_taxa_dict.keys()}
    df['Homolog Tree'] = df.index.map(homolog_tree_dict)

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

    if not new_data and args.ht_include:
        df = update_homolog_tree(df)
    
    df.to_csv(f'{output_fold}/{out_filename}', sep='\t')


def stats_gene(df):
    '''
    Write basic summary for all genes

    :param df: pd.DataFrame with information about all genes and organisms
    :type df: pd.DataFrame
    '''
    # res.write(f'Gene,Total,total[%],SGT\n')
    taxa_count = len(df)
    df[df != 0] = 1
    df = df.sum().to_frame()
    df = df.rename(columns={0: 'Number of Taxa'})
    df[f'Percent of Total Taxa (out of {taxa_count})'] = round((df['Number of Taxa'] / taxa_count) * 100, 2)
    df = df.rename_axis('Gene Name')
    df = df.sort_values(by=['Number of Taxa'], ascending=False)
    df['Homolog Tree'] = ['yes'] * len(df)
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
                          database in downstream homolog tree construction.
                          """))
    optional.add_argument('--ht_include', type=str, metavar='to_include.txt',
                          help=textwrap.dedent("""\
                          Path to text file containing Unique IDs from the database to include in single homolog trees.
                          """))

    in_help = 'Path to fisher.py output directory'
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, in_help=in_help, out_dir=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    args.input_metadata = os.path.abspath(config['PATHS']['input_file'])

    output_fold = os.path.basename(os.path.normpath(args.input)) + '/informant_stats'
    if os.path.isdir(output_fold):
        sys.exit(f'Error: {output_fold} folder already exists.')
    else:
        os.mkdir(output_fold)

    db_taxa_dict = tools.parse_metadata(os.path.join(dfo, 'phylofisher.db'))
    in_taxa_dict = tools.parse_metadata(args.input_metadata, input_meta=True)
    all_taxa_dict = {**db_taxa_dict, **in_taxa_dict}

    matrix, routes = make_table(args.input)

    res = table_with_routes(matrix, routes)

    stats_orgs(matrix)
    stats_orgs(matrix, new_data=True)
    stats_gene(matrix)
