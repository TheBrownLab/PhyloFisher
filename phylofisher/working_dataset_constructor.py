#!/usr/bin/env python
import configparser
import os
from pathlib import Path
from Bio import SeqIO
from phylofisher import help_formatter
from phylofisher.db_map import database, Genes, Sequences


def parse_genes(gene_file):
    '''
    Parse csv table for genes selection

    :param gene_file: path to the gene stats file
    :type gene_file: str
    :return: set of genes to exclude
    :rtype: set
    '''
    to_exlude = set()
    with open(gene_file) as lines:
        next(lines)
        for line in lines:
            gene, _, _, sgt = line.split('\t')
            if sgt.strip().lower() != 'yes':
                to_exlude.add(gene)
    return to_exlude


def parse_orgs(org_file, new_data=False):
    '''
    Parse TSV table for organisms selection

    :param org_file: path to the organism stats file
    :type org_file: str
    :param new_data: Flag indicating if the data is new, defaults to False
    :type new_data: bool, optional
    :return: set of organisms to exclude, set of organisms
    :rtype: tuple(set, set) or set
    '''
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


def fasta_filtr(file, o_to_ex, paralogs=False):
    '''
    Filter fasta sequences with genes. Excludes organisms which should be 
    excluded

    :param file: path to the fasta file
    :type file: str
    :param o_to_ex: set of organisms to exclude
    :type o_to_ex: set
    :param paralogs: Flag indicating if paralogs should be included, defaults to False
    :type paralogs: bool, optional
    '''
    with open(str(Path(args.output, file)), 'w') as res:
        for record in SeqIO.parse(str(Path(args.input, file)), 'fasta'):
            if record.name.split('_')[0] not in o_to_ex:
                res.write(f'>{record.name}\n{record.seq}\n')
        if paralogs:
            # only with paralog selection option
            db_query = Sequences.select(Sequences.header, Sequences.sequence, Sequences.id).join(Genes).where((Genes.name == gene_name) & (Sequences.is_paralog == True))
            
            for q in db_query:
                res.write(f'>{q.header}..p{q.id}\n{q.sequence}\n')


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

    # Connect to database
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    database.init(os.path.join(dfo, 'phylofisher.db'))
    database.connect()

    if args.input[-1] == '/':
        args.input = args.input[:-1]
    os.mkdir(args.output)
    main()

    # Close database connection
    database.close()
