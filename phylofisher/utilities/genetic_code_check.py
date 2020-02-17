#!/usr/bin/env python

import os
import re
import subprocess
from collections import defaultdict
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import AlignIO
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

dataset = '../final_dataset.fas'


def blastDB_tblastn(genome):
    subprocess.call(f'makeblastdb -in {genome} -dbtype nucl -parse_seqids',
                    shell=True)
    subprocess.call(f'tblastn -query {dataset} -db {genome} -num_threads 3 \
                    -outfmt 5 -out output.xml -max_target_seqs 1', shell=True)


def read_transcriptome(transcriptome):
    result_dict = {}
    for sequence in SeqIO.parse(transcriptome, 'fasta'):
        result_dict[sequence.name] = sequence.seq
    return result_dict


def get_sequence_number(file, query):
    sequences = SeqIO.parse(file, 'fasta')
    number = 0
    res = None
    for sequence in sequences:
        if query not in sequence.description:
            number += 1
        elif query in sequence.description:
            res = number
    return res


def alignment_parsing(alignment_file, position):
    # print(alignment_file)
    alignment = AlignIO.read(alignment_file, 'fasta')
    position_dict = defaultdict(list)
    aligment_position = 1
    true_position = 1
    for i in range(alignment.get_alignment_length()):
        if alignment[position, i] != '-':
            position_dict[true_position].append(aligment_position)
            position_dict[true_position].append(alignment[:, i])
            true_position += 1
        aligment_position += 1
    return position_dict


def stops_against_alignment(blast_output, genome):
    genome_dict = {}
    for sequence in SeqIO.parse(genome, 'fasta'):
        genome_dict[sequence.name] = sequence.seq

    blastout = open(blast_output)
    blast_records = NCBIXML.parse(blastout)

    AA = 'RHKDESTNQCGPAVILMFYW-'

    result_dict = defaultdict(list)

    TGA = open('TGA', 'w')
    TAA = open('TAA', 'w')
    TAG = open('TAG', 'w')

    for record in blast_records:
        try:
            gene_name = record.query.split('_')[-1]
            alignment_name = '../alignments/' + gene_name + '.fasta.aln'
            seq_num = get_sequence_number(alignment_name, record.query.split('@')[0])
            positional_dict = alignment_parsing(alignment_name, seq_num)
            if len(record.alignments) > 0:
                for hsp in record.alignments[0].hsps:
                    if '*' in hsp.sbjct[5:-6]:
                        hit_name = record.alignments[0].hit_id
                        query_start = hsp.query_start
                        hit_start = hsp.sbjct_start
                        hit_end = hsp.sbjct_end
                        frame = hsp.frame[1]
                        hit_seq = str(hsp.sbjct)
                        query_seq = str(hsp.query)
                        middline = str(hsp.match)
                        stop_positions = [s.start() for s in re.finditer('\*', hit_seq[5:-6])]
                        stop_positions = [i + 5 for i in stop_positions]
                        for position in stop_positions:
                            if "-" not in query_seq[position - 3:position + 4]:
                                if int(middline[position - 3:position + 4].count(' ')) < 3:
                                    AA_count = {}
                                    real_position = int(query_start + position
                                                        - int(query_seq[:position].count('-')))

                                    if frame in [1, 2, 3]:
                                        real_genome_position = int(
                                            hit_start + (3 * (position - int(hit_seq[:position].count('-')))) - 1)
                                        codon = (
                                        genome_dict[hit_name][real_genome_position:real_genome_position + 3]).upper()
                                    else:
                                        sequence = genome_dict[hit_name]
                                        real_genome_position = int(len(sequence) - hit_end + (
                                                    3 * (position - int(hit_seq[:position].count('-')))))
                                        revcom = sequence.reverse_complement()
                                        codon = (revcom[real_genome_position:real_genome_position + 3]).upper()
                                    col = str(positional_dict[real_position][1])
                                    for amino_acid in AA:
                                        AA_count[amino_acid] = col.count(amino_acid)
                                    y = pd.Series(AA_count)
                                    if y.max() > len(col) * 0.6:
                                        result_dict[codon].append(y.idxmax())
                                        if codon == "TGA":
                                            TGA.write(record.query + '\n')
                                        elif codon == "TAG":
                                            TAG.write(record.query + '\n')
                                        elif codon == "TAA":
                                            TAA.write(record.query + '\n')
        except IndexError:
            print('error in: ', record.query.split('_')[-1])

    TGA.close()
    TAA.close()
    TAG.close()

    return result_dict


def grande_finale(genome):
    basename = genome.split('.')[0]
    os.mkdir(basename)
    subprocess.call(f'mv {genome} {basename}/', shell=True)
    os.chdir(basename)
    blastDB_tblastn(genome)
    res_list_dict = stops_against_alignment('output.xml', genome)

    with PdfPages('AA_stop_counted.pdf') as pdf:
        AA_ = 'RHKDESTNQCGPAVILMFYW'

        for stop, res_list in res_list_dict.items():
            result = {}
            for i in AA_:
                result[i] = res_list.count(i)
                res = pd.Series(result)
            final_result = res.sort_values(ascending=False)
            final_result.plot.bar()
            plt.title(stop)
            plt.xticks(rotation='horizontal')
            pdf.savefig(bbox_inches='tight')
            plt.close()


for file in os.listdir():
    if '.fasta' in file:
        print(file)
        grande_finale(file)
        os.chdir('../')
