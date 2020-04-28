#!/usr/bin/env python
import configparser
import os
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from phylofisher import help_formatter


def check_metadata():
    """
    Checks to make sure taxa provided in to_collapse.csv are actually in the meta data
    """
    with open(args.input, 'r') as infile:
        to_collapse = set()
        for taxon in [line.strip().split(',')[1:] for line in infile]:
            to_collapse.update(taxon)

    taxon_dict = dict.fromkeys(to_collapse, False)

    with open(metadata, 'r') as infile:
        for line in infile:
            line = line.strip()

            short_name = line.split('\t')[0]
            if short_name in to_collapse:
                taxon_dict[short_name] = True

    for taxon in taxon_dict.keys():
        if taxon_dict[taxon] is False:
            print(f'{taxon} is not in the dataset.')

    if False in taxon_dict.values():
        sys.exit('Taxa not collapsed. Please fix errors above.')


def collapse():
    ortho_dir = f'{dfo}/orthologs'
    gene_files = [os.path.join(ortho_dir, gene) for gene in os.listdir(ortho_dir) if gene.endswith('.fas')]
    with open(args.input, 'r') as infile:
        collapse_dict = dict()
        for line in infile:
            collapse_dict[line.strip().split(',')[0]] = line.strip().split(',')[1:]
    empty_record = SeqRecord(Seq('', IUPAC.protein),
                             id='')

    for gene_file in gene_files:
        record_dict = {k: empty_record for k in collapse_dict.keys()}
        records = []
        with open(gene_file, 'r') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                for key in collapse_dict.keys():
                    if record.id in collapse_dict[key]:
                        print('yes')
                        if len(record.seq) > len(record_dict[key].seq):
                            record.id = key
                            record.description = ''
                            record_dict[key] = record
                    else:
                        records.append(record)

        for record in record_dict.values():
            records.append(record)

        with open(gene_file, 'w') as outfile:
            SeqIO.write(records, outfile, 'fasta')


if __name__ == '__main__':
    description = 'Collapses dataset entries into a single taxon'
    parser, optional, required = help_formatter.initialize_argparse(name='taxon_collapser.py',
                                                                    desc=description,
                                                                    usage='taxon_collapser.py [OPTIONS] -i <taxa>')

    in_help = ('Path to CSV formatted file containing taxa to be collapsed.\n'
               'Example Row: "CollapsedName,Taxon1,Taxon2,...,TaxonN')
    args = help_formatter.get_args(parser, optional, required, in_help=in_help, out_dir=False, pre_suf=False)

    # Congig Parser
    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
    metadata = str(os.path.join(dfo, 'metadata.tsv'))

    check_metadata()
    collapse()
