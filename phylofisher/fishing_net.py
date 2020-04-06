#!/usr/bin/env python
import configparser
import os
import textwrap
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
            gene, _, _, sgt = line.split(',')
            sgt = sgt.strip()
            if sgt.lower() != 'yes':
                to_exlude.add(gene)
    return to_exlude


def parse_orgs(org_file):
    """Parse csv table for organisms selection
    input: csv file with organisms
    return: set of organisms to exclude, set of organisms
    for which we want to select paralogs"""
    to_exlude = set()
    paralogs = set()
    with open(org_file) as lines:
        header = next(lines)
        if header.count(',') == 10:
            # only with paralog selection option
            sgt_idx = 9
        else:
            sgt_idx = 6
        for line in lines:
            org = line.split(',')[0]
            sgt = line.split(',')[sgt_idx].strip()
            if sgt.lower() != 'yes':
                to_exlude.add(org)
            if header.count(',') == 10:
                # only with paralog selection option
                if line.split(',')[sgt_idx + 1].strip().lower() == "yes":
                    if org not in to_exlude:
                        paralogs.add(org)
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
    gene_file = str(Path(args.input, 'informant_stats/genes_stats.csv'))
    orgs_file = str(Path(args.input, 'informant_stats/orgs_stats.csv'))
    # genes to exclude
    g_to_ex = parse_genes(gene_file)
    # organisms to exlude
    o_to_ex, paralogs = parse_orgs(orgs_file)
    # only genes which are not in g_to_ex
    filtered_genes = []
    for file in [file for file in os.listdir(args.input) if file.endswith(f"{args.suffix}.fas")]:
        if file.split('.')[0] not in g_to_ex:
            filtered_genes.append(file)

    for file in filtered_genes:
        if args.orthologs is True:
            # Option to only include orthologs
            fasta_filtr(os.path.basename(file), o_to_ex)
        else:
            # Default: Inclueds orthologs and paralogs
            fasta_filtr(os.path.basename(file), o_to_ex, paralogs)


if __name__ == '__main__':
    parser, optional, required = help_formatter.initialize_argparse(name='fishing_net.py',
                                                                    desc='Script for filtering organisms [and|or] genes',
                                                                    usage='fishing_net.py [OPTIONS] '
                                                                          '-i <input> ')

    # Add Arguments Specific to this script
    # Optional
    optional.add_argument('--orthologs', action='store_true',
                          help=textwrap.dedent("""\
                              Only for ortholog selection. Without information
                              about used path."""))

    in_help = "Path to output directory of fisher.py containing the output of informant.py"
    args = help_formatter.get_args(parser, optional, required, in_help=in_help)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['dataset_folder']).resolve())

    if args.input[-1] == '/':
        args.input = args.input[:-1]
    os.mkdir(args.output)
    main()
