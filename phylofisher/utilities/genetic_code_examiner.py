#!/usr/bin/env python
import configparser
import os
import shutil
import subprocess
import textwrap
from collections import defaultdict
from glob import glob
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Data import CodonTable
from matplotlib.backends.backend_pdf import PdfPages

from phylofisher import help_formatter

plt.style.use('ggplot')


def prepare_alignments(threads):
    """
    Prepare multiple sequence alignments from orthologs from the dataset.
    """
    alignmnets_folder = str(Path(dfo, 'alignments'))
    if not os.path.isdir(alignmnets_folder):
        os.mkdir(alignmnets_folder)
    for file in glob(str(Path(dfo, f'orthologs/*.fas'))):
        cmd = f'mafft --auto --reorder --anysymbol --thread {threads} \
            {file} > {file.replace("orthologs", "alignments")}.aln'
        print(cmd)
        subprocess.call(cmd, shell=True)


def parse_query(alignment_file, queries_list):
    """
    Get query sequence (according to input queries) for a given gene. 
    """
    gene = alignment_file.split('/')[-1].split('.')[0]
    ali_dict = {}
    n = 0
    for record in SeqIO.parse(alignment_file, 'fasta'):
        seq = str(record.seq).replace('-', '')
        ali_dict[record.name] = f'>{record.name}@{gene}@{n}\n{seq}\n'
        n += 1
    for query in queries_list:
        if query in ali_dict:
            return ali_dict[query]
    return None


def collect_queries(queries_list):
    """
    Collect query sequences which will be used for blast searches.
    """
    with open('queries.fasta', 'w') as res:
        for file in glob(str(Path(dfo, f'alignments/*.aln'))):
            parsed_query = parse_query(file, queries_list)
            if parsed_query:
                res.write(parsed_query)


def blastDB_tblastn(trans, threads, evalue):
    """
    Make blastdb and perform tblastn with queries from the dataset against given transcriptome.
    """
    subprocess.call(f'makeblastdb -in {trans} -dbtype nucl',
                    shell=True)
    subprocess.call(f'tblastn -query queries.fasta -db {trans} -num_threads {threads} \
                    -evalue {evalue} -outfmt 5 -out blast_out.xml -max_target_seqs 5 ', shell=True)


def read_transcriptome(transcriptome):
    """
    Parse transcriptome as a dictionary.
    """
    result_dict = {}
    for sequence in SeqIO.parse(transcriptome, 'fasta'):
        result_dict[sequence.name] = sequence.seq
    return result_dict


def alignment_parsing(alignment_file, query_n):
    """
    Read multiple sequence aligment and keep information about alignment
    and true positions for all positions for a givne query sequence.
    """
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


def collect_counts(transcriptome, conservation_lvl):
    """
    For every blast record (every gene) check all good quality positions fron alignments
    (avoid first 5 and last 5 amino acids from the alignment; check only
    positions with less than 3 mismatches in their close proximity [-3:+4]) 
    Query position opposite to good qulity position from transcriptome is investigated in
    multiple sequence alignment (prepared only from orthologs from the dataset). In a case,
    when some amino acid is conserved in more then conservation_lvl (70% by default) organims,
    information about conserved amino acid is stored for a corresponding codon in investigated
    transcriptome.
    return: dict of codons and corresponding conserved amino acids dict[UAG] = [Q,Q,Q,E,Q,Q] 
    """

    transcriptome_dict = {}
    for sequence in SeqIO.parse(transcriptome, 'fasta'):
        transcriptome_dict[sequence.name] = sequence.seq

    blastout = open('blast_out.xml')
    blast_records = NCBIXML.parse(blastout)

    AA = 'RHKDESTNQCGPAVILMFYW-'

    result_dict = defaultdict(list)

    for record in blast_records:
        gene_name = record.query.split('@')[1]
        alignment_name = str(Path(dfo, f'alignments/{gene_name}.fas.aln'))
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
                positions = [i + 5 for i in range(len(hit_seq[5:-5]))]
                for position in positions:
                    if "-" not in query_seq[position - 3:position + 4]:
                        if int(middline[position - 3:position + 4].count(' ')) < 3:
                            AA_count = {}
                            real_position = int(query_start + position
                                                - int(query_seq[:position].count('-')))

                            if frame in [1, 2, 3]:
                                real_trans_position = int(
                                    hit_start + (3 * (position - int(hit_seq[:position].count('-')))) - 1)
                                codon = (
                                    transcriptome_dict[hit_name][real_trans_position:real_trans_position + 3]).upper()
                            else:
                                sequence = transcriptome_dict[hit_name]
                                real_trans_position = int(len(sequence) - hit_end + (
                                        3 * (position - int(hit_seq[:position].count('-')))))
                                revcom = sequence.reverse_complement()
                                codon = (revcom[real_trans_position:real_trans_position + 3]).upper()
                            col = str(positional_dict[real_position][1])
                            for amino_acid in AA:
                                AA_count[amino_acid] = col.count(amino_acid)
                            y = pd.Series(AA_count)
                            if y.max() > len(col) * conservation_lvl:
                                result_dict[str(codon)].append(y.idxmax())

    return result_dict


def genecode_plot(res_list_dict, all_codons, transcriptome):
    """
    Prepare bar plot for all codons which seem to differ from the standard  genetic code.
    """
    std_code = CodonTable.unambiguous_dna_by_id[1].forward_table
    std_code['TAG'] = "*"
    std_code['TAA'] = "*"
    std_code['TGA'] = "*"

    with PdfPages(transcriptome + '_genecode.pdf') as pdf:
        AA_ = 'RHKDESTNQCGPAVILMFYW'

        for codon, res_list in res_list_dict.items():
            most_freq = (0, None)
            result = {}
            for i in AA_:
                count = res_list.count(i)
                if count > most_freq[0]:
                    most_freq = (count, i)
                result[i] = res_list.count(i)
                res = pd.Series(result)
            if all_codons:
                if (codon in std_code) and (most_freq[1] != std_code[codon]):
                    print(
                        f"Suspicious codon:{codon.replace('T', 'U')} doens't seem to match The standart genetic code: {std_code[codon]}")
                final_result = res.sort_values(ascending=False)
                final_result.plot.bar()
                plt.title(f"{codon.replace('T', 'U')} (The Standard Code: {std_code[codon]})")
                plt.xticks(rotation='horizontal')
                pdf.savefig(bbox_inches='tight')
                plt.close()
            else:
                if (codon in std_code) and (most_freq[1] != std_code[codon]):
                    print(
                        f"Suspicious codon:{codon.replace('T', 'U')} doens't seem to match The standart genetic code: {std_code[codon]}")
                    final_result = res.sort_values(ascending=False)
                    final_result.plot.bar()
                    plt.title(f"{codon.replace('T', 'U')} (The Standard Code: {std_code[codon]})")
                    plt.xticks(rotation='horizontal')
                    plt.ylabel('Conserved Positions')
                    plt.xlabel('IUPAC Amino Acid Code')
                    pdf.savefig(bbox_inches='tight')
                    plt.close()


def check_alignments():
    """
    Check it aligments folder exist and is not empty.
    """
    alignmnets_folder = str(Path(dfo, 'alignments'))
    if os.path.isdir(alignmnets_folder) == False:
        print("ERROR: For the first run of genetic_code_examiner.py you have to use --prepare_alignments option.")
        return False

    if len(os.listdir(alignmnets_folder)) == 0:
        print("ERROR: Aligments folder is empty. Try to run this script again using --prepare_alignments option.")
        return False
    
    return True

        

def main():
    if args.prepare_alignments:
        prepare_alignments(args.threads)
    if check_alignments(): # execute the block only if alignments are ok
        transcriptome = args.input
        folder = f"{transcriptome.split('.')[0]}_genecode"
        os.mkdir(folder)
        subprocess.call(f'cp {transcriptome} {folder}/', shell=True)
        os.chdir(folder)
        queries_list = args.queries.split(',')
        collect_queries(queries_list)
        print(f'{"=" * 20}\nPerforming blast searches ...\n{"=" * 20}\n')
        blastDB_tblastn(transcriptome, args.threads, args.blast_evalue)
        print(f'{"=" * 20}\nCollecting information from conserved positions ...\n{"=" * 20}\n')
        stats = collect_counts(transcriptome, args.conserved)
        genecode_plot(stats, args.all_codons, transcriptome)
        os.chdir('..')
        if not args.keep_tmp:
            subprocess.call(f'mv {folder}/{transcriptome}_genecode.pdf .', shell=True)
            shutil.rmtree(folder)


if __name__ == '__main__':
    description = 'Script for fast genetic code analysis.'
    usage = 'genetic_code_examiner.py [OPTIONS] -i, --input file.fas -q, --queries Allomacr,Mantplas,...'
    parser, optional, required = help_formatter.initialize_argparse(name='genetic_code_examiner.py',
                                                                    desc=description,
                                                                    usage=usage)

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
    optional.add_argument('-c', '--conserved', type=float, metavar='0-1', default=0.7,
                          help=textwrap.dedent("""\
                          Conservation level. 0-1.
                          Default: 0.7 (70%%)
                          """))
    optional.add_argument('-e', '--blast_evalue', type=str, metavar='1e-X', default='1e-30',
                          help=textwrap.dedent("""\
                          Evalue for blast searches. Default: 1e-30.
                          Default: 1e-30
                          """))
    optional.add_argument('-a', '--all_codons', action='store_true', default=False,
                          help=textwrap.dedent("""\
                          Plot conserved positions for all codons.
                          """))
    optional.add_argument('--keep_tmp', action='store_true',
                          help=textwrap.dedent("""\
                          Keep temporary files
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['database_folder']).resolve())

    main()
