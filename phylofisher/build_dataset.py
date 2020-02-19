#!/usr/bin/env python
import os
import subprocess
import shutil
import argparse
from glob import glob
from collections import defaultdict, Counter
from Bio import SeqIO


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
                res.write(f'>{gene}@{record.name}\n{record.seq}\n')

def diamond():
    """
    Uses for_diamond.fasta to diamond orthomcl.diamonddb

    :return: None
    """

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
        full_name = sline[0] # ADK2@Allomacr
        gene, _ = full_name.split('@') # gene: ADK2
        og = sline[1].split('|')[2].strip() # parse og name OG5_128398
        gene_ogs[gene].append(og) # gene_ogs['ADK2'].append(og)
    return gene_ogs


def get_og_file(threshold):
    prepare_diamond_input() # prepares orthologs for diamond
    diamond() # starts diamond
    gene_filtered_ogs = defaultdict(list)
    gene_ogs = parse_diamond_output()
    for gene, ogs in gene_ogs.items():
        total = len(ogs)
        counted = dict(Counter(ogs))
        for og, count in counted.items():
            if count/total >= threshold:
                gene_filtered_ogs[gene].append(og)
    with open('orthomcl/gene_og', 'w') as res:
        for gene, ogs in gene_filtered_ogs.items():
            res.write(f'{gene}\t{",".join(ogs)}\n')
    os.remove('diamond.res')
    os.remove('for_diamond.fasta')


def datasetdb():
    # TODO fix me
    os.mkdir('datasetdb')
    os.chdir('orthologs')
    subprocess.run('cat *.fas > datasetdb.fasta', shell=True)
    shutil.move('datasetdb.fasta', '../datasetdb')
    os.chdir('../datasetdb')
    dmd_db = 'diamond makedb --in datasetdb.fasta -d prot'
    subprocess.run(dmd_db, shell=True)
    os.remove('datasetdb.fasta')
    os.chdir('..')


def make_profiles(threads):
    """
    Prepares hmm profiles from all fasta files with orthologs
    from orthologs folder.

    :param threads: number of threads
    :return: None
    """


    os.mkdir('profiles')
    os.chdir('orthologs')

    # prepares alignments for all fasta files with orthologs
    aln = f'for i in $(ls *.fas); do mafft --auto' \
        f' --thread {threads} --reorder $i > $(basename $i .fas).aln; done'
    subprocess.run(aln, shell=True)

    # prepares hmm profiles from alignmets made in prev. step
    hmm = 'for i in $(ls *.aln); do hmmbuild $(basename $i .aln).hmm $i; done'
    subprocess.run(hmm, shell=True)
    subprocess.run('mv *.hmm ../profiles', shell=True)
    subprocess.run('rm *.aln', shell=True) # delete alignments
    os.chdir('..')


def main(threads, make_og_file, threshold):
    """
    :param threads: number of threads
    :param make_og_file: Boolean
    :param threshold: float 0-1
    :return: None
    """

    # checks if paralogs folder exists and creates
    # it if not
    if os.path.isdir('paralogs') is False:
        os.mkdir('paralogs')

    datasetdb() # creates datasetdb

    make_profiles(threads)

    if make_og_file:
        get_og_file(threshold)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script for dataset construction.',
                                     usage="build_dataset.py [OPTIONS]")
    parser.add_argument('-t', '--threads', type=int,
                        help='Number of threads, default:1', default=1)
    parser.add_argument('-m', '--make_og_file', action='store_true',
                        help='Make file with information about gene: ogs')
    parser.add_argument('-o', '--og_threshold', type=float,
                        help='Threshold 0-1 for OG. Use only with --make_og_file option.',
                        default=0.1)
    args = parser.parse_args()

    main(args.threads, args.make_og_file, args.og_threshold)
