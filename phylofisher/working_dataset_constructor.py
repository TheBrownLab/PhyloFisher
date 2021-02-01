#!/usr/bin/env python
import configparser
import os
from pathlib import Path

from Bio import SeqIO

from phylofisher import help_formatter


def parse_genes(gene_file):
    """Select gene to exclude: where sgt column is not 'yes' 
    input: csv file with genes
    return: set of genes which should be exluded"""
    to_exlude = set()
    with open(gene_file) as lines:
        next(lines)
        for line in lines:
            gene, _, _, sgt = line.split('\t')
            if sgt.strip().lower() != 'yes':
                to_exlude.add(gene)
    return to_exlude


def parse_orgs(org_file, new_data=False):
    """Parse csv table for organisms selection
    input: csv file with organisms
    return: set of organisms to exclude, set of organisms
    for which we want to select paralogs"""
    to_exlude = set()
    paralogs = set()
    with open(org_file) as lines:
        lines.readline()
        for line in lines:
            if new_data:
                org, _, _, _, _, _, _, _, _, sgt = line.split('\t')
            else:
                org, _, _, _, _, sgt, para = line.split('\t')

            if sgt.strip().lower() != 'yes':
                to_exlude.add(org)

            if not new_data:
                if para.strip().lower() == 'yes':
                    if org not in to_exlude:
                        paralogs.add(org)

        if new_data:
            return to_exlude
        else:
            return to_exlude, paralogs


def fasta_filtr(file, o_to_ex, paralogs=None):
    """Filter fasta sequences with genes. Exludes organisms which should be 
    excluded (o_to_ex)
    input: fasta file with genes, set with orgs to exlude, paralogs=True/None
    result: None"""
    with open(str(Path(args.output, file)), 'w') as res:
        for record in SeqIO.parse(str(Path(args.input, file)), 'fasta'):
            if record.name.split('_')[0] not in o_to_ex:
                res.write(f'>{record.name}\n{record.seq}\n')
        if paralogs:
            # only with paralog selection option
            para_file = str(Path(dfo, f'paralogs/{file.split(".")[0]}_paralogs.fas'))
            if os.path.isfile(para_file):
                for record in SeqIO.parse(para_file, 'fasta'):
                    if record.name.split('.')[0] in paralogs:
                        res.write(f'>{record.name}\n{record.seq}\n')


def main():
    gene_file = str(Path(args.input, 'informant_stats/gene_stats.tsv'))
    new_orgs_file = str(Path(args.input, 'informant_stats/new_taxa_stats.tsv'))
    db_orgs_file = str(Path(args.input, 'informant_stats/db_taxa_stats.tsv'))
    # genes to exclude
    g_to_ex = parse_genes(gene_file)
    # organisms to exlude
    new_o_to_ex = parse_orgs(new_orgs_file, new_data=True)
    db_o_to_ex, paralogs = parse_orgs(db_orgs_file)
    o_to_ex = new_o_to_ex | db_o_to_ex

    # only genes which are not in g_to_ex
    filtered_genes = []
    for file in [file for file in os.listdir(args.input) if file.endswith(f".fas")]:
        if file.split('.')[0] not in g_to_ex:
            filtered_genes.append(file)

    for file in filtered_genes:
        fasta_filtr(os.path.basename(file), o_to_ex, paralogs)


if __name__ == '__main__':
    parser, optional, required = help_formatter.initialize_argparse(name='working_dataset_constructor.py',
                                                                    desc='Script for filtering organisms [and|or] genes',
                                                                    usage='working_dataset_constructor.py [OPTIONS] '
                                                                          '-i <input> ')

    in_help = "Path to output directory of fisher.py containing the output of informant.py"
    args = help_formatter.get_args(parser, optional, required, in_help=in_help, pre_suf=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())

    if args.input[-1] == '/':
        args.input = args.input[:-1]
    os.mkdir(args.output)
    main()
