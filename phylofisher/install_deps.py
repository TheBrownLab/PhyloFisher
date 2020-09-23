#!/usr/bin/env python
import hashlib
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import textwrap
import time
import urllib.request

from phylofisher import help_formatter


def reporthook(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    sys.stdout.write("\r%d%%, %d MB, %d KB/s, %d seconds elapsed" %
                    (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()


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
    urllib.request.urlretrieve(url, f'{fisher_dir}/{fname}', reporthook)

    return fname


def extract(fname):
    tar = tarfile.open(fname, "r:gz")
    tar.extractall()
    tar.close()
    os.remove(fname)


def get_trimal():
    os.chdir(fisher_dir)
    if is_in_path('trimal') is False:
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
        os.chdir(f'{fisher_dir}/standard-RAxML-8.2.12/')
        bash('make -f Makefile.AVX.PTHREADS.gcc')
        # Symlink executables to user bin
        src = f'{fisher_dir}/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-AVX'
        des = f'{user_bin}/raxmlHPC-PTHREADS-AVX2'
        shutil.copy(src, des)


def get_hmmer():
    os.chdir(fisher_dir)
    if is_in_path('hmmsearch') is False:
        # Download and extract hmmer
        os.chdir(f'{fisher_dir}/hmmer-3.3/')
        bash('./configure --prefix=$HOME/.local && make install')


def get_diamond():
    os.chdir(fisher_dir)
    if is_in_path('diamond') is False:
        # Download and extract Diamond
        shutil.copy(f'{fisher_dir}/diamond', f'{user_bin}/diamond')
        os.remove(f'{fisher_dir}/diamond')


def get_fasttree():
    os.chdir(fisher_dir)
    if is_in_path('fasttree') is False:
        # Download and extract Diamond
        bash('gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm')
        shutil.copy(f'{fisher_dir}/FastTree', f'{user_bin}/fasttree')
        os.remove(f'{fisher_dir}/FastTree.c')


def get_blast():
    os.chdir(fisher_dir)
    if is_in_path('blastn') is False:
        if platform.system() == 'Darwin':
            extract(f'{fisher_dir}/ncbi-blast-2.10.1-macosx.tar.gz')
        else:
            extract(f'{fisher_dir}/ncbi-blast-2.10.1-linux.tar.gz')

        src = f'{fisher_dir}/ncbi-blast-2.10.1+/bin/*'
        bash(f'cp {src} {user_bin}')


def get_cd_hit():
    os.chdir(fisher_dir)
    if is_in_path('cd-hit') is False:
        os.chdir(f'{fisher_dir}/cdhit-4.8.1')
        if platform.system() == 'Darwin':
            bash(f'make clean && make CC={args.gxx}')
        else:
            bash('make clean && make')
        exes = ['cd-hit', 'cd-hit-est', 'cd-hit-2d', 'cd-hit-est-2d', 'cd-hit-div', 'cd-hit-454']
        for exe in exes:
            src = f'{fisher_dir}/cdhit-4.8.1/{exe}'
            des = f'{user_bin}/{exe}'
            shutil.copy(src, des)


def get_mafft():
    os.chdir(fisher_dir)
    if is_in_path('mafft') is False:
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
        os.chdir(f'{fisher_dir}/Divvier-1.01')
        bash('make clean && make')
        src = f'{fisher_dir}/Divvier-1.01/divvier'
        des = f'{user_bin}/divvier'
        shutil.copy(src, des)


def get_prequal():
    os.chdir(fisher_dir)
    if is_in_path('prequal') is False:
        os.chdir(f'{fisher_dir}/prequal-1.02')
        bash('make clean && make')
        src = f'{fisher_dir}/prequal-1.02/prequal'
        des = f'{user_bin}/prequal'
        shutil.copy(src, des)


def get_bmge():
    os.chdir(fisher_dir)
    if is_in_path('BMGE') is False:
        os.chdir(f'{fisher_dir}/BMGE-1.12')
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
        des = f'{user_bin}/bmge'
        shutil.copy(src, des)


def get_dist_est():
    os.chdir(fisher_dir)
    if is_in_path('dist_est') is False:
        os.chdir(f'{fisher_dir}/dist_est-1.1')
        bash('make clean && make')
        src = f'{fisher_dir}/dist_est-1.1/dist_est'
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

    # home = os.path.expanduser('~')
    home = os.getcwd()

    user_bin = f'{home}/.local/bin'
    fisher_dir = f'{home}/.local/lib/phylofisher'

    # Check if necessary directy structure is already present.
    # If not creates it
    if os.path.isdir(f'{home}/.local') is False:
        os.mkdir(f'{home}/.local')
    if os.path.isdir(f'{home}/.local/lib') is False:
        os.mkdir(f'{home}/.local/lib')
    if os.path.isdir(user_bin) is False:
        os.mkdir(user_bin)
    if os.path.isdir(fisher_dir) is False:
        os.mkdir(fisher_dir)

    os.chdir(fisher_dir)

    # Get tarball with deps from msstate library
    print('Downloading dependencies:')
    url = 'https://ir.library.msstate.edu/bitstream/handle/11668/19731/Tice_etal.PhyloFisher.archives.tar.gz'
    file_name = 'Tice_etal.PhyloFisher.archives.tar.gz'

    if not os.path.isfile(file_name):
        download(url)

    # Correct original md5 goes here
    original_md5 = '062c84f27388304094d35946f5616a79'

    # Open,close, read file and calculate MD5 on its contents
    with open(file_name, 'rb') as file_to_check:
        # read contents of the file
        data = file_to_check.read()
        # pipe contents of the file through
        md5_returned = hashlib.md5(data).hexdigest()

    # Finally compare original MD5 with freshly calculated
    if original_md5 != md5_returned:
        print('download corrupted')
    else:
        extract(file_name)
        files = os.listdir(f'{fisher_dir}/archives')
        for file in files:
            shutil.move(f'{fisher_dir}/archives/{file}', f'{fisher_dir}/{file}')
            if 'blast' not in file:
                extract(f'{fisher_dir}/{file}')

        shutil.rmtree(f'{fisher_dir}/archives')



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
    get_dist_est()
