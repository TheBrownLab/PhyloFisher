#!/usr/bin/env python

import os
import re
import subprocess
import configparser
import textwrap
from pathlib import Path
from collections import defaultdict
from glob import glob
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import AlignIO
from Bio.Data import CodonTable
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from phylofisher import help_formatter
plt.style.use('ggplot')


def prepare_alignments():
    pass


def parse_query(alignment_file, queries_list):
    gene = alignment_file.split('/')[-1].split('.')[0]
    ali_dict = {}
    n = 0
    for record in SeqIO.parse(alignment_file, 'fasta'):
        seq = str(record.seq).replace('-','')
        ali_dict[record.name] = f'>{record.name}@{gene}@{n}\n{seq}\n'
        n += 1
    for query in queries_list:
        if query in ali_dict:
            return ali_dict[query]
    return None


def collect_queries(queries_list):
    with open('queries.fasta', 'w') as res:
        for file in glob(str(Path(dfo, f'alignments/*.aln'))):
            parsed_query = parse_query(file, queries_list)
            if parsed_query:
                res.write(parsed_query)


def blastDB_tblastn(genome, threads):
    subprocess.call(f'makeblastdb -in {genome} -dbtype nucl',
                    shell=True)
    subprocess.call(f'tblastn -query queries.fasta -db {genome} -num_threads {threads} \
                    -evalue 1e-30 -outfmt 5 -out blast_out.xml -max_target_seqs 5 \
                    -gapopen 11 -gapextend 2', shell=True)


def read_transcriptome(transcriptome):
    result_dict = {}
    for sequence in SeqIO.parse(transcriptome, 'fasta'):
        result_dict[sequence.name] = sequence.seq
    return result_dict


def alignment_parsing(alignment_file, query_n):
    alignment = AlignIO.read(alignment_file, 'fasta')
    position_dict = defaultdict(list)
    aligment_position = 1
    true_position = 1
    for i in range(alignment.get_alignment_length()):
        if alignment[query_n, i] != '-':
            position_dict[true_position].append(aligment_position)
            position_dict[true_position].append(alignment[:, i])
            true_position += 1
        aligment_position += 1
    return position_dict


def stops_against_alignment(transcriptome):
    transcriptome_dict = {}
    for sequence in SeqIO.parse(transcriptome, 'fasta'):
        transcriptome_dict[sequence.name] = sequence.seq

    blastout = open('blast_out.xml')
    blast_records = NCBIXML.parse(blastout)

    AA = 'RHKDESTNQCGPAVILMFYW-'

    result_dict = defaultdict(list)

    for record in blast_records:
        gene_name = record.query.split('@')[1]
        alignment_name = str(Path(dfo, f'alignments/{gene_name}.aln'))
        seq_num = int(record.query.split('@')[-1])
        positional_dict = alignment_parsing(alignment_name, seq_num)
        if len(record.alignments) > 0:
            for hsp in record.alignments[0].hsps:
                hit_name = record.alignments[0].hit_def.split()[0]
                query_start = hsp.query_start
                hit_start = hsp.sbjct_start
                hit_end = hsp.sbjct_end
                frame = hsp.frame[1]
                hit_seq = str(hsp.sbjct)
                query_seq = str(hsp.query)
                middline = str(hsp.match)
                positions = [i + 5 for i in range(len(hit_seq[5:-6]))]
                for position in positions:
                    if "-" not in query_seq[position - 3:position + 4]:
                        if int(middline[position - 3:position + 4].count(' ')) < 3:
                            AA_count = {}
                            real_position = int(query_start + position
                                                - int(query_seq[:position].count('-')))

                            if frame in [1, 2, 3]:
                                real_genome_position = int(
                                    hit_start + (3 * (position - int(hit_seq[:position].count('-')))) - 1)
                                codon = (
                                transcriptome_dict[hit_name][real_genome_position:real_genome_position + 3]).upper()
                            else:
                                sequence = transcriptome_dict[hit_name]
                                real_genome_position = int(len(sequence) - hit_end + (
                                            3 * (position - int(hit_seq[:position].count('-')))))
                                revcom = sequence.reverse_complement()
                                codon = (revcom[real_genome_position:real_genome_position + 3]).upper()
                            col = str(positional_dict[real_position][1])
                            for amino_acid in AA:
                                AA_count[amino_acid] = col.count(amino_acid)
                            y = pd.Series(AA_count)
                            if y.max() > len(col) * 0.7:
                                result_dict[str(codon)].append(y.idxmax())

    return result_dict


def genecode_plot(res_list_dict):
    std_code = CodonTable.unambiguous_dna_by_id[1].forward_table
    std_code['TAG'] = "*"
    std_code['TAA'] = "*"
    std_code['TGA'] = "*"

    with PdfPages('result.pdf') as pdf:
        AA_ = 'RHKDESTNQCGPAVILMFYW*'

        for codon, res_list in res_list_dict.items():
            most_freq = (0, None)
            result = {}
            for i in AA_:
                count = res_list.count(i)
                if count > most_freq[0]:
                    most_freq = (count, i)
                result[i] = res_list.count(i)
                res = pd.Series(result)
            if (codon in std_code) and (most_freq[1] != std_code[codon]):
                final_result = res.sort_values(ascending=False)
                final_result.plot.bar()
                plt.title(f"{codon.replace('T','U')} (The Standard Code: {std_code[codon]})")
                plt.xticks(rotation='horizontal')
                pdf.savefig(bbox_inches='tight')
                plt.close()


def main(args):
    transcriptome = args.input
    folder = f"{transcriptome.split('.')[0]}_genecode"
    os.mkdir(folder)
    subprocess.call(f'cp {transcriptome} {folder}/', shell=True)
    os.chdir(folder)
    queries_list = args.queries.split(',')
    collect_queries(queries_list)
    blastDB_tblastn(transcriptome, args.threads)
    res_list_dict = stops_against_alignment(transcriptome)
    genecode_plot(res_list_dict)

if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['dataset_folder']).resolve())

    description = 'Script for fast genetic code analysis.'
    parser, optional, required = help_formatter.initialize_argparse(name='genetic_code.py',
                                                                    desc=description,
                                                                    usage='genetic_code.py [OPTIONS]',
                                                                    dataset=False,
                                                                    input_meta=False)

    # Required Arguments
    required.add_argument('-i', '--input', required=True, type=str, metavar='file.fas',
                          help=textwrap.dedent("""\
                              Fasta file with nucleotide sequences.
                              """))
    required.add_argument('-q', '--queries', required=True, type=str, metavar='Allomacr,Mantplas,...',
                          help=textwrap.dedent("""\
                              Comma separated short names of organisms which should be used as queries.
                              """))
    # Optional Arguments
    optional.add_argument('-t', '--threads', type=int, metavar='N', default=1,
                        help=textwrap.dedent("""\
                    Number of threads, where N is an integer.
                    Default: 1
                    """))

    optional.add_argument('--prepare_alignments', action='store_true',
                    help=textwrap.dedent("""\
                Prepare alignments for genetic code analysis.
                """))

    optional.add_argument('-c', '--conserved', type=int, metavar='N',
                help=textwrap.dedent("""\
            Conservation level. 0-1.
            """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp=False)
    main(args)