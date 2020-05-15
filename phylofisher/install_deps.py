#!/usr/bin/env python

import os
import platform
import sys
import shutil
import tarfile
import textwrap
import urllib.request
import subprocess
from zipfile import ZipFile

from phylofisher import help_formatter


def bash(cmd):
    subprocess.run(cmd, shell=True)


def is_in_path(cmd):
    """
    Checks to see if command provided is in PATH already. If it is return True, if not returns False
    """
    if shutil.which(cmd) is None:
        in_path = False
    else:
        in_path = True

    return in_path


def download(url):
    fname = url.split('/')[-1]
    urllib.request.urlretrieve(url, f'{fisher_dir}/{fname}')

    return fname


def extract(fname):
    if fname.endswith("tar.gz") or fname.endswith('.tgz'):
        tar = tarfile.open(fname, "r:gz")
        tar.extractall()
        tar.close()
        os.remove(fname)
    elif fname.endswith("tar"):
        tar = tarfile.open(fname, "r:")
        tar.extractall()
        tar.close()
        os.remove(fname)
    elif fname.endswith("zip"):
        zipped = ZipFile(fname, 'r')
        zipped.extractall()
        zipped.close()
        os.remove(fname)


def get_trimal():
    os.chdir(fisher_dir)
    if is_in_path('trimal') is False:
        url = 'http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz'
        extract(download(url))
        os.chdir(f'{fisher_dir}/trimAl/source/')
        bash('make')
        os.chdir(fisher_dir)
        # Symlink executables to user bin
        src = f'{fisher_dir}/trimAl/source/trimal'
        des = f'{user_bin}/trimal'
        shutil.copy(src, des)
        src = f'{fisher_dir}/trimAl/source/readal'
        des = f'{user_bin}/readal'
        shutil.copy(src, des)


def get_raxml():
    os.chdir(fisher_dir)
    if is_in_path('raxmlHPC-PTHREADS-AVX2') is False:
        # Download and extract RAxML
        url = 'https://github.com/stamatak/standard-RAxML/archive/master.zip'
        extract(download(url))
        os.chdir(f'{fisher_dir}/standard-RAxML-master/')
        bash('make -f Makefile.AVX.PTHREADS.gcc')
        # Symlink executables to user bin
        src = f'{fisher_dir}/standard-RAxML-master/raxmlHPC-PTHREADS-AVX'
        des = f'{user_bin}/raxmlHPC-PTHREADS-AVX2'
        shutil.copy(src, des)


def get_hmmer():
    os.chdir(fisher_dir)
    if is_in_path('hmmsearch') is False:
        # Download and extract hmmer
        url = 'http://eddylab.org/software/hmmer/hmmer.tar.gz'
        extract(download(url))

        os.chdir(f'{fisher_dir}/hmmer-3.3/')
        bash('./configure --prefix=$HOME && make install')


def get_diamond():
    os.chdir(fisher_dir)
    if is_in_path('diamond') is False:
        # Download and extract Diamond
        if platform.system() == 'Darwin':
            url = 'http://github.com/bbuchfink/diamond/releases/download/v0.9.31/diamond-linux64.tar.gz'
        else:
            url = 'http://github.com/bbuchfink/diamond/releases/download/v0.9.31/diamond-linux64.tar.gz'
        extract(download(url))
        shutil.copy(f'{fisher_dir}/diamond', f'{user_bin}/diamond')


def get_fasttree():
    os.chdir(fisher_dir)
    if is_in_path('fasttree') is False:
        # Download and extract Diamond
        url = 'http://www.microbesonline.org/fasttree/FastTree.c'
        download(url)
        bash('gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm')
        shutil.copy(f'{fisher_dir}/FastTree', f'{user_bin}/fasttree')


def get_blast():
    os.chdir(fisher_dir)
    if is_in_path('blastn') is False:
        linux_url = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-linux.tar.gz'
        mac_url = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.0+-x64-macosx.tar.gz'
        if platform.system() == 'Darwin':
            extract(download(mac_url))
        else:
            extract(download(linux_url))

        src = f'{fisher_dir}/ncbi-blast-2.10.0+/bin/*'
        bash(f'cp {src} {user_bin}')


def get_cd_hit():
    os.chdir(fisher_dir)
    if is_in_path('cd-hit') is False:
        url = 'https://github.com/weizhongli/cdhit/archive/master.zip'
        extract(download(url))

        os.chdir(f'{fisher_dir}/cdhit-master')
        if platform.system() == 'Darwin':
            bash(f'make clean && make CC={args.gxx}')
        else:
            bash('make clean && make')
        exes = ['cd-hit', 'cd-hit-est', 'cd-hit-2d', 'cd-hit-est-2d', 'cd-hit-div', 'cd-hit-454']
        for exe in exes:
            src = f'{fisher_dir}/cdhit-master/{exe}'
            des = f'{user_bin}/{exe}'
            shutil.copy(src, des)


def get_mafft():
    os.chdir(fisher_dir)
    if is_in_path('mafft') is False:
        # Download and extract mafft
        url = 'https://mafft.cbrc.jp/alignment/software/mafft-7.453-without-extensions-src.tgz'
        extract(download(url))

        os.chdir(f'{fisher_dir}/mafft-7.453-without-extensions/core')
        with open('Makefile', 'r') as infile, open('tmp', 'w') as outfile:
            for line in infile:
                if line.strip() == 'PREFIX = /usr/local':
                    line = f'PREFIX = {fisher_dir}\n'
                elif line.strip() == 'BINDIR = $(PREFIX)/bin':
                    line = f'BINDIR = {user_bin}\n'
                outfile.write(f'{line}')
        print('done')
        shutil.move('tmp', 'Makefile')
        bash('make clean && make && make install')
        os.chdir(fisher_dir)


def get_divvier():
    os.chdir(fisher_dir)
    if is_in_path('divvier') is False:
        url = 'https://github.com/simonwhelan/Divvier/archive/master.zip'
        extract(download(url))
        os.chdir(f'{fisher_dir}/Divvier-master')
        bash('make clean && make')
        src = f'{fisher_dir}/Divvier-master/divvier'
        des = f'{user_bin}/divvier'
        shutil.copy(src, des)


def get_prequal():
    os.chdir(fisher_dir)
    if is_in_path('prequal') is False:
        url = 'https://github.com/simonwhelan/Prequal/archive/master.zip'
        extract(download(url))
        os.chdir(f'{fisher_dir}/prequal-master')
        bash('make clean && make')
        src = f'{fisher_dir}/prequal-master/prequal'
        des = f'{user_bin}/prequal'
        shutil.copy(src, des)


def get_bmge():
    os.chdir(fisher_dir)
    if is_in_path('BMGE') is False:
        url = 'ftp://ftp.pasteur.fr/pub/gensoft/projects/BMGE/BMGE-1.12.tar.gz'
        bash(f'wget {url}')
        extract('BMGE-1.12.tar.gz')
        os.chdir('BMGE-1.12')
        with open('stub.sh', 'w') as outfile:
            stub = ('#!/bin/sh\n'
                    'MYSELF=`which "$0" 2>/dev/null`\n'
                    '[ $? -gt 0 -a -f "$0" ] && MYSELF="./$0"\n'
                    'java=java\n'
                    'if test -n "$JAVA_HOME"; then\n'
                    '    java="$JAVA_HOME/bin/java"\n'
                    'fi\n'
                    'java_args=-Xmx1g\n'
                    'exec "$java" $java_args -jar $MYSELF "$@"\n'
                    'exit 1\n')
            outfile.write(stub)

        bash('cat stub.sh BMGE.jar > BMGE && chmod +x BMGE ')
        src = f'{fisher_dir}/BMGE-1.12/BMGE'
        des = f'{user_bin}/BMGE'
        shutil.copy(src, des)


def get_dist_est():
    os.chdir(fisher_dir)
    if is_in_path('dist_est') is False:
        url = 'https://github.com/DavidZihala/PhyloFisher/raw/master/lib/dist_estv1.1.tar.gz'
        extract(download(url))
        os.chdir(f'{fisher_dir}/dist_estv1.1')
        bash('make clean && make')
        src = f'{fisher_dir}/dist_estv1.1/dist_est'
        des = f'{user_bin}/dist_est'
        shutil.copy(src, des)


if __name__ == '__main__':
    description = 'Downloads, compiles, and installs trimal, RAxML, hmmer, diamond, FastTree, BLAST+, cd-hit' \
                  'mafft, Divvier, prequal, BMGE, and dist_est'
    parser, optional, required = help_formatter.initialize_argparse(name='install_deps.py',
                                                                    desc=description,
                                                                    usage='install_deps.py [OPTIONS]')

    optional.add_argument('--gxx', type=str, metavar='<path/to/g++>', default=None,
                          help=textwrap.dedent("""\
                            FOR OSX ONLY: Path to g++ compiler.
                            If not provided on OSX, cd-hit installation will be skipped 
                            """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    home = os.path.expanduser('~')
    user_bin = f'{home}/bin'
    fisher_dir = f'{home}/PhyloFisher_lib'

    if os.path.isdir(user_bin) is False:
        os.mkdir(user_bin)
    if os.path.isdir(fisher_dir) is False:
        os.mkdir(fisher_dir)

    os.chdir(fisher_dir)
    get_trimal()
    get_raxml()
    get_hmmer()
    get_diamond()
    get_fasttree()
    get_blast()
    if platform.system() == 'Darwin' and args.gxx is None:
        print('Your operating system is OSX, and a g++ compiler was not provided.\n'
              'cd-hit will NOT be installed')
    else:
        get_cd_hit()
    get_mafft()
    get_divvier()
    get_prequal()
    get_bmge()
    # get_dist_est()
