#!/usr/bin/env python
import os
import subprocess
import textwrap
from collections import defaultdict
from glob import glob
from multiprocessing import Pool
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofisher import help_formatter

# SNAKEFILE_PATH = f'{os.path.dirname(os.path.realpath(__file__))}/nucl_matrix_constructor.smk'

def bash(cmd):
    '''
    Run bash command

    :param cmd: bash command
    :type cmd: str
    '''
    subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)


def get_genes():
    '''
    Get input gene files

    :return: input gene files
    :rtype: list
    '''
    genes = [file.split(".")[0] for file in os.listdir(args.input)]

    return genes


def parse_input_tsv(input_tsv):
    '''
    Parse input tsv file

    :param input_tsv: input tsv file
    :type input_tsv: str
    :return: dictionary of unique IDs and paths to nucleotide files
    :rtype: dict
    '''
    ret = defaultdict(dict)
    with open(input_tsv, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            ret[line[0]] = line[1]

    return ret


def make_blast_db(fasta_dict):
    '''
    Make blast database

    :param fasta_dict: dictionary of unique IDs and paths to nucleotide files
    :type fasta_dict: dict
    '''
    for k in fasta_dict:
        bash(f'makeblastdb -in {fasta_dict[k]} -dbtype nucl')


def tblastn(fasta_dict, genes):
    '''
    Run tblastn

    :param fasta_dict: dictionary of unique IDs and paths to nucleotide files
    :type fasta_dict: dict
    :param genes: list of genes
    :type genes: list
    '''
    os.makedirs(os.path.dirname(f'{args.output}/tblastn'), exist_ok=True)

    for gene in genes:
        with open(f'{args.input}/{gene}.fas', 'r') as infile:
            for record in SeqIO.parse(infile, 'fasta'):
                bash(f'echo -e ">{record.description}\n{record.seq}" | tblastn -db {fasta_dict[record.description]} -out {args.output}/tblastn/{gene}.{record.description}.tsv -outfmt "6 sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads {args.threads} -max_target_seqs 1 -max_hsps 1 -evalue 1e-10')


def parse_blast(fasta_dict):
    os.makedirs(os.path.dirname(f'{args.output}/nucl_seqs'), exist_ok=True)
    
    seq_dict = {}    
    for file in os.listdir(f'{args.output}/tblastn'):
        taxon = os.path.basename(file).split('.')[1]
        gene = os.path.basename(file).split('.')[0]
        with open(f'{args.output}/tblastn/{gene}.{taxon}.tsv', 'r') as infile, open(fasta_dict[taxon], 'r') as nucl_file:
            for line in infile:
                sseqid, _, _, _, _, _, _, sstart, send, _, _ = line.strip().split('\t')

            for record in SeqIO.parse(nucl_file, 'fasta'):
                if record.id == sseqid:
                    record.seq = record.seq[int(sstart):int(send)]
                    record.description = taxon
                    if gene not in seq_dict:
                        seq_dict[gene] = [record]
                    else:
                        seq_dict[gene].append(record)
    
    for gene in seq_dict:
        with open(f'{args.output}/nucl_seqs/{gene}.fas', 'w') as outfile:
            SeqIO.write(seq_dict[gene], outfile, 'fasta')


if __name__ == '__main__':
    description = 'To get nucleotides trim align and concatenate orthologs into a nucleotide super-matrix:'
    parser, optional, required = help_formatter.initialize_argparse(name='nucl_matrix_constructor.py',
                                                                    desc=description,
                                                                    usage='nucl_matrix_constructor.py [OPTIONS] -i path/to/input/')

    # Required Arguments
    required.add_argument('-it', '--input_tsv', metavar='<path>', type=str, required=True,
                          help=textwrap.dedent("""\
                          Path to input tsv file which contains unique IDs and paths to nucleotide files.
                          """))
    # Optional Arguments
    optional.add_argument('-of', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Desired format of the output matrix.
                          Options: fasta, phylip (names truncated at 10 characters), 
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Format of the input files.
                          Options: fasta, phylip (names truncated at 10 characters), 
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))
    optional.add_argument('-t', '--threads', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Desired number of threads to be utilized.
                          Default: 1
                          """))
    optional.add_argument('--clean_up', action='store_true', default=False,
                          help=textwrap.dedent("""\
                          Clean up large intermediate files.
                          """))

    # Changes help descriptions from the default input and output help descriptions
    in_help = 'Path to prep_final_dataset_<M.D.Y>'
    args = help_formatter.get_args(parser, optional, required, in_help=in_help)

    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)


    fasta_dict = parse_input_tsv(args.input_tsv)
    # make_blast_db(fasta_dict)
    # tblastn(fasta_dict, get_genes())
    parse_blast(fasta_dict)