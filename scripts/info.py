import glob
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import numpy as np
import argparse


def taxonomy_dict(metadata, multi_input):
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
            if record.name.count('_') == 3:
                name, _, path, _ = record.name.split('_')
                paths[gene_name][name] = path
            else:
                name = record.name
            names.add(name)
    return sorted(list(names)), paths


def get_gene_column(gene, names):
    gene_name = gene.split('/')[-1].split('.')[0]
    column = pd.Series(np.zeros(len(names)), index=names, name=gene_name)
    for record in SeqIO.parse(gene, "fasta"):
        if record.name.count('_') == 3:
            org = record.name.split('_')[0]
        else:
            org = record.name
        column[org] += 1
    return column.astype(int)


def make_table(folder):
    genes = glob.glob(f'{folder}/*.fas')
    names, paths = collect_names(genes)
    columns = []
    for gene in genes:
        columns.append(get_gene_column(gene, names))
    df = pd.DataFrame(columns)
    df = df.transpose()
    df = df.reindex(sorted(df.columns), axis=1)
    df.to_csv("alvert_test.csv")
    return df, paths


def table_with_paths(df, paths):
    tax_list = []
    for org in df.index:
        try:
            org_tax = t_dict[org]
        except KeyError:
            org_tax = "Unknown"
        tax_list.append(org_tax)
    tax_col = pd.Series(tax_list, df.index)
    df.insert(loc=0, column='Taxonomy', value=tax_col)
    for gene in df.columns:
        df[gene] = df[gene].apply(str)
        for org in df[gene].index:
            if org in paths[gene]:
                df.at[org, gene] = f'{df[gene][org]}_{paths[gene][org]}'
    print(df)
    return df


# TODO org taxonomy genes genes_missing genes_with_SBH genes_with_BBH gene_with_hmm
def stats_orgs(table):
    rows = []
    for org in table.index:
        genes_tot = len(table.columns) - 1
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
        for val in res.loc[org].values[1:]:
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
    df.to_csv("alvert2_test.csv")
    return df


def stats_gene(table):
    columns = table.iloc[:, 1:]
    tab_len = len(table)
    with open("gene_stats", 'w') as res:
        for i in columns:
            orgs = tab_len - table[i].value_counts()['0']
            orgs_perc = (orgs/tab_len) * 100
            res.write(f"{i},{orgs},{orgs_perc}\n")



t_dict, fnames = taxonomy_dict('/home/david/MsPhylo/data/metadata.tsv', '/home/david/MsPhylo/test/test_input2')
tab, paths = make_table("/home/david/MsPhylo/test/Apr28/fasta")
res = table_with_paths(tab, paths)
res.to_csv('presence_path.csv')
stats_orgs(res)
stats_gene(res)