#!/usr/bin/env python
import os
import subprocess
import shutil
import argparse


def datasetdb():
    os.mkdir('datasetdb')
    os.chdir('orthologs')
    subprocess.run('cat *.fas > datasetdb.fasta', shell=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    shutil.move('datasetdb.fasta', '../datasetdb')
    os.chdir('../datasetdb')
    dmd_db = 'diamond makedb --in datasetdb.fasta -d prot'
    subprocess.run(dmd_db, shell=True)
    os.remove('datasetdb.fasta')
    os.chdir('..')


def make_profiles(threads):
    os.mkdir('profiles')
    os.chdir('orthologs')
    aln = f'for i in $(ls *.fas); do mafft --auto' \
        f' --thread {threads} --reorder $i > $(basename $i .fas).aln; done'
    subprocess.run(aln, shell=True)

    hmm = 'for i in $(ls *.aln); do hmmbuild $(basename $i .aln).hmm $i; done'
    subprocess.run(hmm, shell=True)
    subprocess.run('mv *.hmm ../profiles', shell=True)
    subprocess.run('rm *.aln', shell=True)


def main(threads):
    if os.path.isdir('paralogs') is False:
        os.mkdir('paralogs')

    datasetdb()
    make_profiles(threads)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script for dataset construction.',
                                     usage="build_dataset.py [OPTIONS]")
    parser.add_argument('-t', '--threads', type=int,
                        help='Number of threads, default:1', default=1)
    args = parser.parse_args()

    main(args.threads)