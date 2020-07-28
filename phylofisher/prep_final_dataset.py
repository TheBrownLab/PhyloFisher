import configparser
import os
import shutil
import textwrap
from pathlib import Path

from Bio import SeqIO

from phylofisher import help_formatter


def parse_ortholog_tsv():
    """

    :return:
    """
    with open('select_orthologs.tsv', 'r') as infile:
        infile.readline()
        genes_to_include = []
        for line in infile:
            line = line.strip()
            taxon, _, include = line.split('\t')
            if include == 'yes':
                genes_to_include.append(taxon)

    return genes_to_include


def subset_orthologs():
    if os.path.isdir(args.output) is False:
        os.makedirs(args.output)

    genes = parse_ortholog_tsv()
    files = [os.path.join(orthologs_dir, x) for x in os.listdir(orthologs_dir) if x.endswith('.fas')]
    for gene in genes:
        for file in files:
            if gene == os.path.basename(file).split('.')[0]:
                src = file
                dest = f'{args.output}/{os.path.basename(file)}'
                shutil.copy(src, dest)


def parse_taxa_tsv():
    """

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


def subset_taxa():
    """

    :return:
    """
    taxa = parse_taxa_tsv()

    files = [file for file in os.listdir(args.output)]
    for file in files:
        with open(f'{args.output}/{file}', 'r') as infile, open(f'{args.output}/tmp', 'w') as outfile:
            records = []
            for record in SeqIO.parse(infile, 'fasta'):
                if record.description in taxa:
                    records.append(record)

            SeqIO.write(records, outfile, 'fasta')

        shutil.move(f'{args.output}/tmp', f'{args.output}/{file}')


if __name__ == '__main__':
    description = 'Subsets taxa and orthologs to be included in super matrix construction'
    parser, optional, required = help_formatter.initialize_argparse(name='prep_final_dataset.py',
                                                                    desc=description,
                                                                    usage='prep_final_dataset.py '
                                                                          '[OPTIONS]')

    # Optional Arguments
    # optional.add_argument('--orthologs', type=str, metavar='orthologs.tsv', default=None,
    #                       help=textwrap.dedent("""\
    #                           Path to orthologs.tsv. Output of select_orthologs.py.
    #                           """))
    # optional.add_argument('--taxa', type=str, metavar='taxa.tsv', default=None,
    #                       help=textwrap.dedent("""\
    #                           Path to taxa.tsv. Output of select_taxa.py
    #                           """))

    args = help_formatter.get_args(parser, optional, required, inp_dir=False, pre_suf=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
    orthologs_dir = f'{dfo}/orthologs'

    if os.path.isfile('select_orthologs.tsv'):
        subset_orthologs()
    else:
        files = [os.path.basename(file) for file in os.listdir(orthologs_dir)]
        for file in files:
            shutil.copy(f'{orthologs_dir}/{file}', f'{args.output}/{file}')

    if os.path.isfile('select_taxa.tsv'):
        subset_taxa()
