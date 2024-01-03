#!/usr/bin/env python
import csv
import os
import sys
import subprocess
import tarfile
import textwrap
from collections import defaultdict
from glob import glob
from multiprocessing import Pool
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofisher import help_formatter

out_dict = {'fasta':          'fas',
            'phylip':         'phy',
            'phylip-relaxed': 'phy',
            'nexus':          'nex'}

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
    Parse metadata

    :param input_tsv: metadata
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

def cp_cds_files(fasta_dict):
    os.makedirs(f'{args.output}/cds', exist_ok=True)

    for k in fasta_dict:
        # Determine if tar.gz and open appropriately
        if fasta_dict[k].endswith('.tar.gz'):
            with tarfile.open(fasta_dict[k], 'r:gz') as tar:
                with tar.extractfile(fasta_dict[k].split('.tar.gz')[0]) as fasta_file:
                    records = []
                    for i, record in enumerate(SeqIO.parse(fasta_file, 'fasta')):
                        record.id = ''
                        record.description = f'{k}_{i}'
                        records.append(record)
            
            SeqIO.write(records, outfile, 'fasta')
            fasta_dict[k] = f'{args.output}/cds/{k}.fas'
        else:  
            with open(fasta_dict[k], 'r') as infile:
                records = []
                for i, record in enumerate(SeqIO.parse(infile, 'fasta')):
                    record.id = ''
                    record.description = f'{k}_{i}'
                    records.append(record)
        
        with open(f'{args.output}/cds/{k}.fas', 'w') as outfile:
            SeqIO.write(records, outfile, 'fasta')
            fasta_dict[k] = f'{args.output}/cds/{k}.fas'
    
    return fasta_dict

def make_blast_db(k):
    '''
    Make blast database

    :param fasta_dict: dictionary of unique IDs and paths to nucleotide files
    :type fasta_dict: dict
    '''
    os.makedirs(f'{args.output}/logs/makeblastdb', exist_ok=True)
    bash(f'makeblastdb -in {fasta_dict[k]} -dbtype nucl &> {args.output}/logs/makeblastdb/{k}.log')


def tblastn(gene):
    '''
    Run tblastn

    :param fasta_dict: dictionary of unique IDs and paths to nucleotide files
    :type fasta_dict: dict
    :param genes: list of genes
    :type genes: list
    '''
    # make output directory and log directory
    os.makedirs(f'{args.output}/tblastn', exist_ok=True)
    os.makedirs(f'{args.output}/logs/tblastn', exist_ok=True)

    # tblastn parameters
    out_fmt = '6 sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    max_target_seqs = 1
    max_hsps = 1
    evalue = '1e-10'

    with open(f'{args.input}/{gene}.fas', 'r') as infile:
        for record in SeqIO.parse(infile, 'fasta'):
            cmd_frags = (
                f'echo -e ">{record.description}\n{record.seq}" |',
                f'tblastn', 
                    f'-db {fasta_dict[record.description]}',
                    f'-out {args.output}/tblastn/{gene}.{record.description}.tsv',
                    f'-outfmt "{out_fmt}"',
                    f'-num_threads {args.threads}',
                    f'-max_target_seqs {max_target_seqs}',
                    f'-max_hsps {max_hsps}',
                    f'-evalue {evalue}',
                    f'&> {args.output}/logs/tblastn/{gene}.{record.description}.log'
            )
            
            bash(' '.join(cmd_frags))

def parse_blast():
    '''
    Get nucleotide sequences from cds files

    :param fasta_dict: dictionary of unique IDs and paths to nucleotide files
    :type fasta_dict: dict
    '''
    os.makedirs(f'{args.output}/nucl_seqs', exist_ok=True)
    
    seq_dict = {}

    for file in os.listdir(f'{args.output}/tblastn'):
        taxon = os.path.basename(file).split('.')[1]
        gene = os.path.basename(file).split('.')[0]
        # open blast results and cds file
        with open(f'{args.output}/tblastn/{gene}.{taxon}.tsv', 'r') as infile, open(fasta_dict[taxon], 'r') as nucl_file:
            # loop through blast results
            for line in infile:
                sseqid, _, _, _, _, _, _, sstart, send, _, _ = line.strip().split('\t')
            
            # loop through cds file
            for record in SeqIO.parse(nucl_file, 'fasta'):
                if sseqid in record.id + record.description:
                    # get nucleotide sequence                    
                    if int(sstart) > int(send):
                        record.seq = record.seq[int(send):int(sstart)]
                        record.seq = record.seq.reverse_complement()
                    else:
                        record.seq = record.seq[int(sstart):int(send)]
                    
                    # Set record id and description to taxon name
                    record.id = taxon
                    record.description = ''
                    
                    # add to dictionary
                    if gene not in seq_dict:
                        seq_dict[gene] = [record]
                    else:
                        seq_dict[gene].append(record)
                    
                    # break loop once found
                    break
    
    # write to fasta file
    for gene in seq_dict:
        with open(f'{args.output}/nucl_seqs/{gene}.fas', 'w') as outfile:
            SeqIO.write(seq_dict[gene], outfile, 'fasta')


def align(gene):
    '''
    _summary_

    :param fasta_dict: _description_
    :type fasta_dict: _type_
    '''
    # make output and log directories 
    os.makedirs(f'{args.output}/mafft', exist_ok=True)
    os.makedirs(f'{args.output}/logs/mafft', exist_ok=True)
    os.makedirs(f'{args.output}/trimal', exist_ok=True)
    os.makedirs(f'{args.output}/logs/trimal', exist_ok=True)
    
    gene = os.path.basename(gene).split('.')[0]
    
    # mafft
    cmd_frags = (
        'mafft --auto',
            f'--thread {args.threads}',
            f'{args.output}/nucl_seqs/{gene}.fas',
            f'> {args.output}/mafft/{gene}.aln',
            f'2> {args.output}/logs/mafft/{gene}.log'
    )

    bash(' '.join(cmd_frags))

    # trimAl
    # parameters
    gt = 0.8
    # Command fragments
    cmd_frags = (
        'trimal',
            f'-in {args.output}/mafft/{gene}.aln',
            f'-out {args.output}/trimal/{gene}.final',
            f'-gt {gt}',
            f'&> {args.output}/logs/trimal/{gene}.log'
    )

    bash(' '.join(cmd_frags))

def get_orgs(input_folder):
    '''

    :param input_folder:
    :return:
    '''
    name_set = set()
    files = sorted(glob(f'{input_folder}/*'))
    for file in files:
        with open(file) as f:
            for record in SeqIO.parse(f, args.in_format):
                fname = record.description
                name = fname.split('_')[0]
                name_set.add(name)
    return sorted(list(name_set))


def concatenate():
    '''
    '''
    with open(f'{args.output}/indices.tsv', 'w') as outfile:
        outfile.write('Gene\tStart\tStop\n')
        files = sorted(glob(f'{args.output}/trimal/*'))

        orgs = get_orgs(args.input)
        occupancy_dict = {k:[] for k in orgs}

        total_len = 0
        res_dict = defaultdict(str)
        for file in files:
            gene = os.path.basename(file).split('.')[0]
            length = 0
            seq_dict = {}
            myformat = 'fasta'
            for record in SeqIO.parse(file, myformat):
                length = len(record.seq)
                seq_dict[record.id.split('_')[0]] = str(record.seq)
            start_len = total_len + 1
            total_len += length
            outfile.write(f'{gene}\t{start_len}\t{total_len}\n')
            for org in orgs:
                if org in seq_dict:
                    res_dict[org] += seq_dict[org]
                    occupancy_dict[org].append(1)
                else:
                    res_dict[org] += ('-' * length)
                    occupancy_dict[org].append(0)
        
        occupancy_df = pd.DataFrame(occupancy_dict, index=[os.path.basename(file).split('.')[0] for file in files]).transpose()
        occupancy_df.to_csv(f'{args.output}/occupancy.tsv', sep='\t')
    
    records = []
    for org, seq in res_dict.items():
        records.append(SeqRecord(Seq(seq),
                                id=org,
                                name='',
                                description=''))
                                
    if args.out_format.lower() in out_dict:
        with open(f'{args.output}/matrix.{args.out_format.lower()}', "w") as handle:
            SeqIO.write(records, handle, args.out_format.lower())
    else:
        sys.exit('Invalid Output Format')

    with open(f'{args.output}/matrix_constructor_stats.tsv', 'w') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Taxon', 'PercentMissingData'])
        missing = []
        for record in SeqIO.parse(f'{args.output}/matrix.{args.out_format.lower()}', args.out_format.lower()):
            missing.append((record.name, (record.seq.count('-') / total_len) * 100))
        for org_missing in sorted(missing, key=lambda x: x[1], reverse=True):
            tsv_writer.writerow(list(org_missing))


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
    args = help_formatter.get_args(parser, optional, required, in_help=in_help, pre_suf=False)

    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)

    os.makedirs(f'{args.output}', exist_ok=True)

    fasta_dict = parse_input_tsv(args.input_tsv)

    fasta_dict = cp_cds_files(fasta_dict)
    
    with Pool(args.threads) as p:
        p.map(make_blast_db, fasta_dict.keys())

    with Pool(args.threads) as p:
        p.map(tblastn, get_genes())
    
    parse_blast()

    with Pool(args.threads) as p:
        p.map(align, os.listdir(f'{args.output}/nucl_seqs'))
    
    concatenate()