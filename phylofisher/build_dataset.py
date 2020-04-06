#!/usr/bin/env python
import csv
import os
import random
import string
import subprocess
import shutil
import argparse
import textwrap
from glob import glob
from tempfile import NamedTemporaryFile
from collections import defaultdict, Counter
from Bio import SeqIO
import re

from phylofisher import help_formatter


def check_ortho_taxa():
    tempfile = NamedTemporaryFile(mode='w', delete=False)
    # Creates set of taxa in orthologs dir
    unique_orgs_orthos = set()
    files = os.listdir('./orthologs')
    for file in files:
        file = os.path.join('./orthologs', file)
        with open(file, 'r') as infile, tempfile:
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    name = line[1:]
                    # Replaces @'s and ..'s in FASTA headers with _'s
                    if '@' in name:
                        name = name.replace('@', '_')
                    if '..' in name:
                        name = name.replace('..', '_')
                    tempfile.write(f'>{name}\n')
                else:
                    tempfile.write(line + '\n')
                unique_orgs_orthos.add(name)
        # Overwrite old file with new file with new headers
        shutil.move(tempfile.name, file)

        return unique_orgs_orthos


def check_meta_taxa():
    # Creates set of taxa in metadata
    tempfile = NamedTemporaryFile(mode='w', delete=False)
    unique_orgs_meta = set()
    with open('metadata.tsv', 'r') as csvfile, tempfile:
        reader = csv.reader(csvfile, delimiter='\t')
        writer = csv.writer(tempfile, delimiter='\t')
        for row in reader:
            if row[0] == 'Name in Dataset':
                continue
            if '@' in row[0]:
                row[0] = row[0].replace('@', '_')
            if '..' in row[0]:
                row[0] = row[0].replace('..', '_')
            unique_orgs_meta.add(row[0])
            writer.writerow(row)
    # Overwrite old file with new file with new headers
    shutil.move(tempfile.name, 'metadata.tsv')

    return unique_orgs_meta


def check_taxa():
    """
    Checks to see if all taxa in metadata is in ortholog dir and vice-versa. If there are differences the script exits
    here and prints the discrepancies. Also replaces @'s and ..'s with _'s if present in the UniqueID

    :return: None
    """
    unique_orgs_orthos = check_ortho_taxa()
    unique_orgs_meta = check_meta_taxa()

    # Exits if there are differences and prints those differences
    if len(unique_orgs_meta - unique_orgs_orthos) > 0 or len(unique_orgs_orthos - unique_orgs_meta) > 0:
        if len(unique_orgs_meta - unique_orgs_orthos) > 0:
            print(f'Taxa in metadata but not in orthologs dir:')
            for org in (unique_orgs_meta - unique_orgs_orthos):
                print(org)
        if len(unique_orgs_orthos - unique_orgs_meta) > 0:
            print(f'Taxa in orthlog dir but not in metadata:')
            for org in (unique_orgs_orthos - unique_orgs_meta):
                print(org)
        os.sys.exit()


def id_generator(size=5, chars=string.digits):
    """"Generate random number with 5 digits."""
    return ''.join(random.choice(chars) for _ in range(size))


def paralog_name(abbrev, keys):
    """"Prepare paralog name (short name + 5 random digits).
    Recursive function.
    example: Homosap..12345
    input: short name of an organism, names of already existing paralogs
    for a given organism
    return: unique paralog name"""
    id_ = id_generator()
    pname = f'{abbrev}..p{id_}'
    if pname not in keys:
        return pname
    else:
        paralog_name(abbrev, keys)


def check_paralogs():
    """
    Checks formatting of headers in paralog files, and appends "..pNNNNN" to header if not already present
    """
    tempfile = NamedTemporaryFile(mode='w', delete=False)
    # Creates set of taxa in orthologs dir
    unique_orgs_orthos = set()
    files = os.listdir('./paralogs')
    id_set = set()

    for file in files:
        file = os.path.join('./paralogs', file)
        with open(file, 'r') as infile, tempfile:
            for line in infile:
                if line.startswith('>'):
                    name = line.strip()[1:]
                    if re.search(r'\.\.p{5}\d', name):
                        id_set.add(name[-5:])
                    else:
                        name = paralog_name(name, id_set)
                        id_set.add(name[-5:])
                    tempfile.write(f'{name}\n')
                else:
                    tempfile.write(line + '\n')


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
    # os.remove('diamond.res')
    # os.remove('for_diamond.fasta')


def concat_gene_files():
    files = glob('*.fas')
    print()
    for file in files:
        with open(file, 'r') as infile, open('datasetdb.fasta', 'w') as outfile:
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    gene = ''.join(os.path.basename(file).split(".")[:-1])
                    line = f'{line}@{gene}'
                outfile.write(f'{line}\n')


def datasetdb():
    # TODO fix me
    os.mkdir('datasetdb')
    os.chdir('orthologs')
    concat_gene_files()
    shutil.move('datasetdb.fasta', '../datasetdb')
    os.chdir('../datasetdb')
    dmd_db = 'diamond makedb --in datasetdb.fasta -d datasetdb'
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
    subprocess.run('rm *.aln', shell=True)  # delete alignments
    os.chdir('..')


def main(threads, no_og_file, threshold):
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
    else:
        check_paralogs()

    datasetdb()  # creates datasetdb

    make_profiles(threads)

    if not no_og_file:
        get_og_file(threshold)


if __name__ == "__main__":
    description = 'Script for dataset construction.'
    parser, optional, required = help_formatter.initialize_argparse(name='build_dataset.py',
                                                                    desc=description,
                                                                    usage='build_dataset.py [OPTIONS]')

    optional.add_argument('-t', '--threads', type=int,
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
                            Threshold 0-1 for OG.
                            Default: 0.1
                            """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    check_taxa()
    main(args.threads, args.no_og_file, args.og_threshold)
