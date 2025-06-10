#!/usr/bin/env python
import csv
import os
import shutil
import subprocess
import sys
import textwrap
from collections import defaultdict, Counter
from glob import glob
import matplotlib.colors as mcolors
from Bio import SeqIO
from phylofisher import help_formatter
from phylofisher.db_map import database, Genes, Sequences


def prepare_diamond_input():
    """
    Prepares seqs for diamond in a way that it connects
    information about gene before sequence name.
    Example: Allomacr from ADK2.fas will be named
    ADK2@Allomacr at the for_diamond.fasta file.

    :return: None
    """
    files = glob('orthologs/*.fas')
    with open("for_diamond.fasta", 'w') as res:
        for file in files:
            gene = file.split('/')[-1].split('.')[0]
            for record in SeqIO.parse(file, 'fasta'):
                res.write(f'>{gene}@{record.description}\n{record.seq}\n')


def diamond():
    """
    Uses for_diamond.fasta to diamond orthomcl.diamonddb

    :return: None
    """
    with open('for_diamond.fasta', 'r') as infile, open('new.fas', 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line.startswith(">") is False:
                line = line.replace('-', '')
            outfile.write(f'{line}\n')

    shutil.move('new.fas', 'for_diamond.fasta')

    db = 'orthomcl/orthomcl.diamonddb.dmnd'
    out = 'diamond.res'
    cmd = (
        f'diamond blastp -e 1e-10 -q for_diamond.fasta --more-sensitive '
        f'--db {db} -o {out} -p {args.threads} --outfmt 6 qseqid stitle evalue')
    subprocess.run(cmd, shell=True)


def parse_diamond_output():
    gene_ogs = defaultdict(list)
    for line in open('diamond.res'):
        sline = line.split("\t")
        full_name = sline[0]  # ADK2@Allomacr
        gene, _ = full_name.split('@')  # gene: ADK2
        og = sline[1].split('|')[2].strip()  # parse og name OG5_128398
        gene_ogs[gene].append(og)  # gene_ogs['ADK2'].append(og)
    return gene_ogs


def genes_in_orthodb():
    """
    Check if gene has an og in orthomcl. If not -> report this gene.
    """
    gene_in_gene_og = set()
    for line in open("orthomcl/gene_og"):
        gene_in_gene_og.add(line.split()[0])

    files = glob('orthologs/*.fas')
    genes = [file.split('/')[-1].split('.')[0] for file in files]
    rerun = False
    for gene in genes:
        if gene not in gene_in_gene_og:
            print(f"{gene} has not hit in orthomcl and thus can not be used in databases."
            "\nPlease delete this gene in orthologs and re-run build_database.py.")
            rerun = True
    if rerun == True:
        sys.exit()


def get_og_file(threshold):
    prepare_diamond_input()  # prepares orthologs for diamondF
    diamond()  # starts diamond
    gene_filtered_ogs = defaultdict(list)
    gene_ogs = parse_diamond_output()
    for gene, ogs in gene_ogs.items():
        total = len(ogs)
        counted = dict(Counter(ogs))
        for og, count in counted.items():
            if count / total >= threshold:
                gene_filtered_ogs[gene].append(og)
    with open('orthomcl/gene_og', 'w') as res:
        for gene, ogs in gene_filtered_ogs.items():
            res.write(f'{gene}\t{",".join(ogs)}\n')
    os.remove('diamond.res')
    os.remove('for_diamond.fasta')
    

def datasetdb(threads):
    '''
    

    :param threads: threads to use for multithreading
    :type threads: int
    '''
    os.mkdir('datasetdb')
    os.chdir('datasetdb')

    with open('datasetdb.fasta', 'w') as outfile:
        db_query = Sequences.select(
            Sequences.header, Sequences.sequence, Sequences.id).where(
                Sequences.is_paralog == False)
        for q in db_query:
            outfile.write(f'>{q.header}\n{q.sequence}\n')
    
    dmd_db = 'diamond makedb --in datasetdb.fasta -d datasetdb'
    subprocess.run(dmd_db, shell=True)

    os.mkdir('../profiles')
    os.chdir('../profiles')


    # prepares alignments for all fasta files with orthologs
    for q in Genes.select(Genes.id, Genes.name):
        gene_id = q.id
        gene = q.name
        alns = ''
        db_query = Sequences.select(
            Sequences.header, Sequences.sequence, Sequences.id).where(
                Sequences.is_paralog == False, Sequences.gene_id == gene_id)
        
        for q in db_query:
            alns += f'>{q.header}\n{q.sequence}\n'

        with open('tmp.fas', 'w') as outfile:
            outfile.write(alns)
        
        aln = f'mafft --globalpair --thread {threads} --reorder tmp.fas> {gene}.fas.aln'
        subprocess.run(aln, shell=True, check=True)
        os.remove('tmp.fas')  # delete tmp file

        # prepares hmm profiles from alignmets made in prev. step
        hmm = f'hmmbuild {gene}.hmm {gene}.fas.aln'
        subprocess.run(hmm, shell=True, check=True)
        os.remove(f'{gene}.fas.aln') # delete alignment file


def generate_tree_colors():

    colors = sorted(mcolors.CSS4_COLORS.keys())
    higher_tax_set = set()
    with open('metadata.tsv', 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            higher_tax_set.add(row[2])
    
    with open('tree_colors.tsv', 'w') as outfile:
        outfile.write('Taxonomy\tColor\n')
        for i, higher_tax in enumerate(sorted(list(higher_tax_set))):
            outfile.write(f'{higher_tax}\t{colors[i]}\n')

def main(threads, no_og_file, threshold):
    """
    :param args:
    :param threads: number of threads
    :param make_og_file: Boolean
    :param threshold: float 0-1
    :return: None
    """

    if os.path.isdir('datasetdb') is True:
        shutil.rmtree('datasetdb')
    if os.path.isdir('profiles') is True:
        shutil.rmtree('profiles')

    datasetdb(threads)  # creates datasetdb

    if not no_og_file:
        get_og_file(threshold)

    genes_in_orthodb()

if __name__ == "__main__":
    description = 'Script for database construction and taxonomic updates. Must be run within path/to/database/'
    parser, optional, required = help_formatter.initialize_argparse(name='build_database.py',
                                                                    desc=description,
                                                                    usage='build_database.py [OPTIONS]')

    optional.add_argument('-t', '--threads', type=int, default=1,
                          help=textwrap.dedent("""\
                          Number of threads.
                          Default:1
                          """))

    optional.add_argument('-n', '--no_og_file', action='store_true',
                          help=textwrap.dedent("""\
                          Do not make Gene OG file
                          """))
    optional.add_argument('-o', '--og_threshold', type=float, default=0.1, metavar='0.X',
                          help=textwrap.dedent("""\
                          0.X proportion of sequences (0-1) of sequences that must hit an 
                          OrthoMCL orthogroup for the group to be assigned.
                          Default: 0.1 (10%%)
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    main(args.threads, args.no_og_file, args.og_threshold)
