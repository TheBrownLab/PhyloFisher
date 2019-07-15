#!/usr/bin/env python
import glob
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import numpy as np
import argparse
from pathlib import Path
import configparser
import os
import sys


def taxonomy_dict(metadata, multi_input=None):
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
    if multi_input:
        for line in open(multi_input):
            metadata_input = line.split('\t')
            tax = metadata_input[2].strip()
            group = metadata_input[3].strip()
            full_name = metadata_input[5].strip()
            tax_g[tax] = group
            full_names[tax] = full_name
    return tax_g, full_names


def collect_names(files):
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
    df = df.reindex(sorted(df.columns), axis=1)
    return df, paths


def table_with_paths(df, paths):
    full_names = []
    tax_list = []
    for org in df.index:
        try:
            org_tax = t_dict[org]
        except KeyError:
            org_tax = "Unknown"
        tax_list.append(org_tax)
        full_names.append(fnames[org])

    tax_col = pd.Series(tax_list, df.index)
    df.insert(loc=0, column='Taxonomy', value=tax_col)

    fname_col = pd.Series(full_names, df.index)
    df.insert(loc=0, column='full_name', value=fname_col)

    df.to_csv(f'{output_fold}/occupancy.csv')

    for gene in df.columns:
        df[gene] = df[gene].apply(str)
        for org in df[gene].index:
            if org in paths[gene]:
                df.at[org, gene] = f'{df[gene][org]}_{paths[gene][org]}'
    return df


def stats_orgs_path(table):
    #TODO paralogs no/available
    rows = []
    for org in table.index:
        genes_tot = len(table.columns) - 2
        try:
            genes = genes_tot - table.loc[org].value_counts()['0']
        except KeyError:
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
        rows.append(pd.Series([fnames[org], t_dict[org], genes, missing, missing_perc, sbh['SBH'], sbh['BBH'],
                               sbh['HMM'], 'yes', 'no'],
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
    return df


def stats_orgs(table):
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
    return df

def stats_gene(table):
    columns = table.iloc[:, 2:]
    tab_len = len(table)
    with open(f'{output_fold}/genes_stats.csv', 'w') as res:
        res.write(f'gene,total,total[%],SGT\n')
        for i in columns:
            orgs = tab_len - table[i].value_counts()['0']
            orgs_perc = (orgs/tab_len) * 100
            res.write(f"{i},{orgs},{orgs_perc.round(2)},yes\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='informant', usage="informant.py -i input_folder [OPTIONS]")
    parser.add_argument('-i', '--input_folder', required=True)
    parser.add_argument('-m', '--metadata', required=True)
    parser.add_argument('-n', '--input_metadata')
    parser.add_argument('-s', '--sufix')
    parser.add_argument('--paralog_selection', action='store_true')
    parser.add_argument('--occupancy_with_paths', action='store_true')
    # parser.add_argument('--orthologs', action='store_true')
    args = parser.parse_args()


    # config = configparser.ConfigParser()
    # config.read('config.ini')
    # dfo = str(Path(args.metadata).resolve())
    # multi_input = os.path.abspath(config['PATHS']['input_file'])

    # if args.orthologs:
    #     args.input_folder = Path(dfo, 'orthologs')

    output_fold = os.path.basename(Path(args.input_folder)) + '_stats'
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