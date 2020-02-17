# !/usr/bin/env python
import glob
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import numpy as np
import argparse
import configparser
from pathlib import Path
import os
import sys
import textwrap
from phylofisher import fisher


def taxonomy_dict(metadata, input_metadata=None):
    """
    Read metadata from dataset and input_metadata.
    input: metadata file, input metadata file (optional)
    return: dictionary with taxon: group; dictionary with
    full names 
    """
    tax_g = {}
    full_names = {}
    for line_ in open(metadata):
        if 'Full Name' not in line_:
            sline = line_.split('\t')
            tax = sline[0].strip()
            group = sline[2].strip()
            full_name = sline[1].strip()
            tax_g[tax] = group
            full_names[tax] = full_name
    if input_metadata:
        for line in open(input_metadata):
            metadata_input = line.split('\t')
            tax = metadata_input[2].strip()
            group = metadata_input[3].strip()
            full_name = metadata_input[6].strip()
            tax_g[tax] = group
            full_names[tax] = full_name
    return tax_g, full_names


def collect_names(files):
    """
    Collect all short names. Collest information about 'path'
    for newly added sequences (BBH, SBH, HMM).
    input: fasta files with genes (from fisher)
    return: set of sorted shortnames; dictionary with genes 
    as a keys and deeper level dictionaries with shortnames of organisms 
    as keys and BBH or SBH or HMM as values
    """
    names = set()
    paths = defaultdict(dict)
    for file in files:
        gene_name = file.split('/')[-1].split('.')[0]
        for record in SeqIO.parse(file, 'fasta'):
            if record.name.count('_') >= 3:
                name = record.name.split('_')[0]
                path = record.name.split('_')[2]
                paths[gene_name][name] = path
            elif '_' in record.name:
                name = record.name.split('_')[0]
            else:
                name = record.name
            names.add(name)
    return sorted(list(names)), paths


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
    'paths' used for sequence selection (BBH, SBH, HMM)
    """
    if args.sufix:
        genes = glob.glob(f'{folder}/*{args.sufix}')
    else:
        genes = glob.glob(f'{folder}/*')
    names, paths = collect_names(genes)
    columns = []
    for gene in genes:
        columns.append(get_gene_column(gene, names))
    df = pd.DataFrame(columns)
    df = df.transpose()
    df = df.reindex(sorted(df.columns), axis=1)  # comment me mf
    return df, paths


def table_with_paths(df, paths):
    """Create occupancy table and add information about paths used for sequence
    selection (BBH, SBH, HMM) to input dataframe
    input: pd.DataFrame with information about all genes (present(>0)/absent(0)) 
    for all organisms, dictionary with paths for all sequences
    output: modifiend input dataframe with information about path
    """
    # TODO: refractor me please
    full_names = []
    tax_list = []
    for org in df.index:
        try:
            org_tax = t_dict[org]
        except KeyError:
            org_tax = "Unknown"  # why mf?
        tax_list.append(org_tax)
        full_names.append(fnames[org])

    tax_col = pd.Series(tax_list, df.index)
    df.insert(loc=0, column='Taxonomy', value=tax_col)

    fname_col = pd.Series(full_names, df.index)
    df.insert(loc=0, column='full_name', value=fname_col)

    df.to_csv(f'{output_fold}/occupancy.csv')

    # Adds paths to df
    for gene in df.columns:
        df[gene] = df[gene].apply(str)
        for org in df[gene].index:
            if org in paths[gene]:
                df.at[org, gene] = f'{df[gene][org]}_{paths[gene][org]}'
    return df


def paralog_orgs():
    """Collect all shortnames for organisms with at least one paralog
    in the dataset.
    input: None
    return: set of short names of organisms with at leat one paralog
    """
    paralogs = set()
    paralog_fold = os.path.dirname(args.metadata)
    for file in glob.glob(f'{paralog_fold}/paralogs/*.fas'):
        for record in SeqIO.parse(file, 'fasta'):
            paralogs.add(record.name.split('.')[0])
    return paralogs


def stats_orgs_path(table):
    """
    Create csv file with basic summary about analyzed dataset with information
     about used 'paths' (BBH, SBH, HMM) and paralogs.
    input: dataframe
    return: None
    """
    paralogs = paralog_orgs()
    rows = []
    for org in table.index:
        genes_tot = len(table.columns) - 2  # because org name is the index
        try:
            # subtracting genes with value == 0
            genes = genes_tot - table.loc[org].value_counts()['0']
        except KeyError:
            # why mf? Probably no '0'
            genes = genes_tot
        missing = genes_tot - genes
        missing_perc = (missing / genes_tot) * 100
        sbh = {}
        sbh['SBH'] = 0
        sbh['BBH'] = 0
        sbh['HMM'] = 0
        for val in res.loc[org].values[2:]:
            if '_' in val:
                path = val.split('_')[1]
                sbh[path] += 1
        para_ava = 'none'
        if org in paralogs:
            para_ava = 'yes'
        rows.append(pd.Series([fnames[org], t_dict[org], genes, missing, missing_perc, sbh['SBH'], sbh['BBH'],
                               sbh['HMM'], 'yes', para_ava],
                              index=["full name", "taxonomy", "#Genes", "#Missing", '%Missing', "#SBH",
                                     "#BBH", "#HMM", "SGT", "paralogs"],
                              name=org))
    df = pd.DataFrame(rows)
    df["#Genes"] = df['#Genes'].astype(int)
    df["#Missing"] = df['#Missing'].astype(int)
    df["%Missing"] = df['%Missing'].round(2)
    df["#SBH"] = df["#SBH"].astype(int)
    df["#BBH"] = df["#BBH"].astype(int)
    df["#HMM"] = df["#HMM"].astype(int)
    df.to_csv(f'{output_fold}/orgs_stats.csv')


def stats_orgs(table):
    """
    Create csv file with basic summary about analyzed dataset without information
    about 'paths' and paralogs.
    input: dataframe
    return: None
    """
    rows = []
    for org in table.index:
        genes_tot = len(table.columns) - 2
        try:
            genes = genes_tot - table.loc[org].value_counts()['0']
        except KeyError:
            genes = genes_tot
        missing = genes_tot - genes
        missing_perc = (missing / genes_tot) * 100
        rows.append(pd.Series([fnames[org], t_dict[org], genes, missing, missing_perc, "yes"],
                              index=["full name", "taxonomy", "#Genes", "#Missing", '%Missing', "SGT"],
                              name=org))
    df = pd.DataFrame(rows)
    df["#Genes"] = df['#Genes'].astype(int)
    df["#Missing"] = df['#Missing'].astype(int)
    df["%Missing"] = df['%Missing'].round(2)
    df.to_csv(f'{output_fold}/orgs_stats.csv')


def stats_gene(table):
    """Write basic summary for all genes to a file
    input: dataframe
    return: None
    """
    columns = table.iloc[:, 2:]
    tab_len = len(table)
    with open(f'{output_fold}/genes_stats.csv', 'w') as res:
        res.write(f'gene,total,total[%],SGT\n')
        for i in columns:
            orgs = tab_len - table[i].value_counts()['0']
            orgs_perc = (orgs / tab_len) * 100
            res.write(f"{i},{orgs},{orgs_perc.round(2)},yes\n")


if __name__ == '__main__':
    formatter = lambda prog: fisher.myHelpFormatter(prog, max_help_position=100)

    parser = argparse.ArgumentParser(prog='informant.py',
                                     description='some description',
                                     usage='informant.py -i input_folder [OPTIONS]',
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent("""\
                                         additional information:
                                            stuff
                                            """))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    # Required Arguments
    required.add_argument('-i', '--input', required=True, type=str, metavar='path/to/input/',
                          help=textwrap.dedent("""\
                          Path to input directory
                          """))
    required.add_argument('-m', '--metadata', required=True, type=str, metavar='meta.tsv',
                          help=textwrap.dedent("""\
                          Path to metadata.tsv
                          """))

    # Optional Arguments
    optional.add_argument('-n', '--input_metadata', type=str, metavar='new',
                          help=textwrap.dedent("""\
                          Meta data from newly added data
                          """))
    optional.add_argument('--paralog_selection', action='store_true',
                          help=textwrap.dedent("""\
                          some description
                              """))
    optional.add_argument('--orthologs', action='store_true',
                          help=textwrap.dedent("""\
                          some description
                          """))
    optional.add_argument('--occupancy_with_paths', action='store_true',
                          help=textwrap.dedent("""\
                          some description
                          """))
    optional.add_argument('-c', '--use_config', action='store_true',
                          help=textwrap.dedent("""\
                          Use config
                          """))
    optional.add_argument('-s', '--suffix', metavar='<suffix>', type=str,
                          help=textwrap.dedent("""\
                          Suffix of input files
                          Default: NONE
                          Example: path/to/input/*.suffix
                          """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                          Show this help message and exit.
                          """))

    parser._action_groups.append(optional)
    args = parser.parse_args()

    if args.use_config:
        config = configparser.ConfigParser()
        config.read('config.ini')
        dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
        args.input_metadata = os.path.abspath(config['PATHS']['input_file'])
        args.metadata = os.path.join(dfo, 'metadata.tsv')
        if args.orthologs:
            args.input_folder = os.path.join(dfo, 'orthologs')

    output_fold = os.path.basename(os.path.normpath(args.input_folder)) + '_stats'
    if os.path.isdir(output_fold):
        sys.exit(f'Error: {output_fold} folder already exists.')
    else:
        os.mkdir(output_fold)

    t_dict, fnames = taxonomy_dict(args.metadata, args.input_metadata)
    tab, paths = make_table(args.input_folder)
    res = table_with_paths(tab, paths)

    if args.occupancy_with_paths:
        res.to_csv(f'{output_fold}/occupancy_with_paths.csv')
    if args.paralog_selection:
        stats_orgs_path(res)
    else:
        stats_orgs(res)
    stats_gene(res)
