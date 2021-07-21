#!/usr/bin/env python
import configparser
import os
import shutil
import sys
import textwrap
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofisher import help_formatter


def parse_ortholog_tsv():
    """
    Parses select_orthologs.tsv (output of select_orthologs.py) to determine which genes are to be included in
    the final dataset.
    :return:
    """
    with open('select_orthologs.tsv', 'r') as infile:
        header = infile.readline()
        genes_to_include = []
        for line in infile:
            line = line.strip()
            if 'In-Group Completeness' in header:
                gene, _, _, _, include = line.split('\t')
            else:
                gene, _, include = line.split('\t')
            if include == 'yes':
                genes_to_include.append(gene)

    return genes_to_include


def subset_orthologs():
    """
    Subsets orthologs from the Database into a new directory to be included in the final dataset
    :return:
    """
    # Creates output dir if it doesn't already exist
    if os.path.isdir(args.output) is False:
        os.makedirs(args.output)

    # Genes to include in the final dataset
    genes = parse_ortholog_tsv()
    # List of paths to ortholog files
    files = [os.path.join(orthologs_dir, x) for x in os.listdir(orthologs_dir) if x.endswith('.fas')]
    # Copies gene files to the output dir
    for gene in genes:
        for file in files:
            if gene == os.path.basename(file).split('.')[0]:
                src = file
                dest = f'{args.output}/{os.path.basename(file)}'
                shutil.copy(src, dest)


def parse_taxa_tsv():
    """
    Parses select_taxa.tsv (output of select_taxa.py) to determine which taxa are to be included in
    the final dataset.
    :param tsv_file:
    :return:
    """
    with open('select_taxa.tsv', 'r') as infile:
        infile.readline()
        taxa_to_include = []
        for line in infile:
            line = line.strip()
            taxon, _, _, _, include = line.split('\t')
            if include == 'yes':
                taxa_to_include.append(taxon)

    return taxa_to_include


def get_chimeras():
    """
    Parses chimeras.tsv if the chimera flag is utilized
    :return: Dict with chimera's Unique ID as the key and a list the Unique IDs of taxa comprising the chimera as the
    value
    """
    chim_dict = dict()
    chim_taxa = []
    with open(args.chimeras, 'r') as infile:
        for line in infile:
            split_line = line.strip().split('\t')
            chim_dict[split_line[0]] = split_line[3:]
            chim_taxa += split_line[3:]
    return chim_dict, chim_taxa


def subset_taxa():
    """
    Subsets taxa from the Database into a new directory to be included in the final dataset
    """
    taxa = parse_taxa_tsv()

    if args.chimeras:
        chimeras, chimeric_taxa = get_chimeras()

    files = [file for file in os.listdir(args.output)]

    for file in files:
        with open(f'{args.output}/{file}', 'r') as infile, open(f'{args.output}/tmp', 'w') as outfile:
            records = []
            if args.chimeras:
                chimera_seqs = {k: '' for k in chimeras.keys()}

            for record in SeqIO.parse(infile, 'fasta'):
                # Only includes taxa marked "yes" in select_taxa.tsv

                if record.description in taxa:
                    records.append(record)

                elif args.chimeras and record.description in chimeric_taxa:
                    for key in chimeras.keys():
                        if record.description in chimeras[key] and len(record.seq) > len(chimera_seqs[key]):
                            chimera_seqs[key] = str(record.seq)


            if args.chimeras:
                for org, seq in chimera_seqs.items():
                    if len(seq) > 0:
                        records.append(SeqRecord(Seq(seq), id=org, name='', description=''))

            SeqIO.write(records, outfile, 'fasta')

        shutil.move(f'{args.output}/tmp', f'{args.output}/{file}')


def check_if_empty():
    """
    Checks to see if there are empty files in the output directory
    :return: Boolian value (True if no empty files)
    """
    for subdir, dirs, files in os.walk(args.output):
        for file in files:
            filepath = subdir + os.sep + file

            if os.path.getsize(filepath) == 0:
                sys.exit('There are empty files in the output. Please check for errors')

    return True


if __name__ == '__main__':
    description = 'Subsets taxa and orthologs to be included in super matrix construction'
    parser, optional, required = help_formatter.initialize_argparse(name='prep_final_dataset.py',
                                                                    desc=description,
                                                                    usage='prep_final_dataset.py '
                                                                          '[OPTIONS]')

    # Optional Arguments
    optional.add_argument('--chimeras', type=str, metavar='chimeras.tsv', default=None,
                          help=textwrap.dedent("""\
                          Path to chimeras.tsv. This will collapses taxa listed in chimera.tsv into a chimera 
                          keeping the longest sequence for each gene.
                          """))

    args = help_formatter.get_args(parser, optional, required, inp_dir=False, pre_suf=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())
    orthologs_dir = f'{dfo}/orthologs'

    if os.path.isdir(args.output) is False:
        os.mkdir(args.output)

    if os.path.isfile('select_orthologs.tsv'):
        subset_orthologs()
    else:
        files = [os.path.basename(file) for file in os.listdir(orthologs_dir)]
        for file in files:
            shutil.copy(f'{orthologs_dir}/{file}', f'{args.output}/{file}')

    if os.path.isfile('select_taxa.tsv'):
        subset_taxa()

    check_if_empty()
