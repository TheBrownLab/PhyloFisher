#!/usr/bin/env python
import csv
import os
import shutil
import subprocess
import textwrap
from functools import partial
from multiprocessing import Pool
from tempfile import NamedTemporaryFile

from Bio import SeqIO

from phylofisher import help_formatter


def bash_command(cmd):
    """Function to run bash commands in a shell"""
    command_run = subprocess.call(cmd, shell=True, executable='/bin/bash')#,
                                #  stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if command_run == 0:
        return True
    else:
        return False


def mkdir_and_cd(dir_name):
    try:
        os.mkdir(dir_name)
    except FileExistsError:
        pass
    os.chdir(dir_name)


def delete_gaps_stars(root):
    """Removes -'s and *'s from alignments"""
    file = f'{args.input}/{root}.fas'
    file_name = f'{root}.aa'
    with open(file_name, 'w') as res:
        for record in SeqIO.parse(file, 'fasta'):
            res.write(f'>{record.name}\n{str(record.seq).replace("-", "").replace("*", "")}\n')


def x_to_dash(file):
    """Replaces X's in alignments with -'s"""
    file_name = f'{os.path.basename(file).split(".")[0]}.pre_bmge'
    with open(f'{file_name}', 'w') as res:
        for record in SeqIO.parse(file, 'fasta'):
            res.write(f'>{record.name}\n{str(record.seq).replace("X", "-")}\n')


def read_full_proteins(core):
    full_prots = {}
    for record in SeqIO.parse(f'{args.output}/prequal/{core}.aa.filtered', 'fasta'):
        full_prots[record.name] = record.seq
    return full_prots


def good_length(trimmed_aln, threshold):
    core = trimmed_aln.split('.')[0]
    full_proteins = read_full_proteins(core)
    original_name = f'{core}.length_filtered'
    length = None
    with open(original_name, 'w') as res:
        for record in SeqIO.parse(trimmed_aln, 'fasta'):
            if length is None:
                length = len(record.seq)
            coverage = len(str(record.seq).replace('-', '').replace('X', '')) / len(record.seq)
            if coverage > threshold:
                res.write(f'>{record.description}\n{full_proteins[record.name]}\n')
            else:
                print('deleted:', record.name, coverage)


def mk_checkpoint_tmp(files):
    if os.path.isfile(f'checkpoint.tmp') is False:
        with open(f'checkpoint.tmp', 'w') as outfile:
            header = 'taxon,prequal,mafft1,divvier1,bmge,mafft2,divvier2,trimal,raxml\n'
            outfile.write(header)
            for file in files:
                outfile.write(f'{file}{9 * ",0"}\n')
    else:
        pass


def get_checkpoints():
    checks = {}
    if os.path.isfile(f'{args.output}/checkpoint.tmp') is True:
        with open(f'checkpoint.tmp', 'r') as infile:
            for line in infile:
                line = line.strip()
                checks[line.split(',')[0]] = line.split(',')[1:]

    return checks


def update_checkpoints(root, prog):
    filename = f'{args.output}/checkpoint.tmp'
    tempfile = NamedTemporaryFile(mode='w', delete=False)

    fields = ['taxon', 'prequal', 'mafft1', 'divvier1', 'bmge', 'mafft2', 'divvier2', 'trimal', 'raxml']

    with open(filename, 'r') as csvfile, tempfile:
        reader = csv.DictReader(csvfile, fieldnames=fields)
        writer = csv.DictWriter(tempfile, fieldnames=fields)
        # Writes header to the output file
        writer.writerow(next(reader))
        for row in reader:
            if row['taxon'] == str(root):
                row[prog] = 1
            row = {'taxon': row['taxon'], 'prequal': row['prequal'], 'mafft1': row['mafft1'],
                   'divvier1': row['divvier1'], 'bmge': row['bmge'], 'mafft2': row['mafft2'],
                   'divvier2': row['divvier2'], 'trimal': row['trimal'], 'raxml': row['raxml']
                   }
            writer.writerow(row)
    shutil.move(tempfile.name, filename)


def prepare_analyses(checks, root):
    if (args.threads / file_count) > 1:
        threads = int(args.threads / file_count)
    else:
        threads = 1
    checks = checks[root]

    for i, check in enumerate(checks):
        # prequal
        if i == 0 and check == '0':
            mkdir_and_cd(f'{args.output}/prequal')
            delete_gaps_stars(root)
            status = bash_command(f'prequal {root}.aa')
            if status:
                update_checkpoints(root, 'prequal')
            else:
                break

        # Length Filtration
        elif i == 1 and check == '0':
            mkdir_and_cd(f'{args.output}/length_filtration')
            # mafft1
            mkdir_and_cd('mafft')
            # MAFFT - length filtration
            status = bash_command(
                f'mafft --thread {threads} --globalpair --maxiterate 1000 --unalignlevel 0.6 '
                f'{args.output}/prequal/{root}.aa.filtered > {root}.aln')
            if status:
                update_checkpoints(root, 'mafft1')
            else:
                break

        elif i == 2 and check == '0':
            mkdir_and_cd(f'{args.output}/length_filtration')
            mkdir_and_cd(f'{args.output}/length_filtration/divvier')
            # Divvier - length filtration
            status = bash_command(f'divvier -mincol 4 -partial {args.output}/length_filtration/mafft/{root}.aln')
            if status:
                shutil.move(f'../mafft/{root}.aln.divvy.fas', f'./{root}.aln.divvy.fas')
                shutil.move(f'../mafft/{root}.aln.PP', f'./{root}.aln.PP')
                update_checkpoints(root, 'divvier1')
            else:
                break

        elif i == 3 and check == '0':
            mkdir_and_cd(f'{args.output}/length_filtration')
            mkdir_and_cd(f'{args.output}/length_filtration/bmge')
            # outputs pre_bmge
            x_to_dash(f'{args.output}/length_filtration/divvier/{root}.aln.divvy.fas')
            # BMGE
            status = bash_command(f'BMGE -t AA -g 0.3 -i {root}.pre_bmge -of {root}.bmge')
            if status:
                good_length(trimmed_aln=f'{root}.bmge', threshold=0.5)
                update_checkpoints(root, 'bmge')
            else:
                break

        # mafft2
        elif i == 4 and check == '0':
            mkdir_and_cd(f'{args.output}/mafft')
            status = bash_command(
                f'mafft --thread {threads} --globalpair --maxiterate 1000 --unalignlevel 0.6 '
                f'{args.output}/length_filtration/bmge/{root}.length_filtered > {root}.aln2')
            if status:
                update_checkpoints(root, 'mafft2')
            else:
                break

        # divvier2
        elif i == 5 and check == '0':
            mkdir_and_cd(f'{args.output}/divvier')
            status = bash_command(f'divvier -mincol 4 -partial {args.output}/mafft/{root}.aln2')
            if status:
                shutil.move(f'../mafft/{root}.aln2.divvy.fas', f'./{root}.aln2.divvy.fas')
                shutil.move(f'../mafft/{root}.aln2.PP', f'./{root}.aln2.PP')
                update_checkpoints(root, 'divvier2')
            else:
                break

        # trimal
        elif i == 6 and check == '0':
            mkdir_and_cd(f'{args.output}/trimal')
            status = bash_command(f'trimal -in {args.output}/divvier/{root}.aln2.divvy.fas -gt 0.01 -out {root}.final')
            if status:
                update_checkpoints(root, 'trimal')
            else:
                break

        # raxml
        elif args.no_tree is False and i == 7 and check == '0':
            mkdir_and_cd(f'{args.output}/RAxML')
            status = bash_command(f'raxmlHPC-PTHREADS-AVX2 -T {threads} -m PROTGAMMALG4XF -f a '
                                  f'-s {args.output}/trimal/{root}.final -n {root}.tre -x 123 -N 100 -p 12345')
            if status:
                update_checkpoints(root, 'raxml')
            else:
                break

        elif not args.no_tree and i == 8:
            # add_length(root, length)
            mkdir_and_cd(f'{args.output}/trees')
            shutil.copy(f'{args.output}/RAxML/RAxML_bipartitions.{root}.tre',
                        f'{args.output}/trees/RAxML_bipartitions.{root}.tre')
            shutil.copy(f'{args.output}/trimal/{root}.final',
                        f'{args.output}/trees/{root}.final')
            shutil.copy(f'{args.output}/length_filtration/bmge/{root}.bmge',
                        f'{args.output}/trees/{root}.trimmed')

    return root, checks


if __name__ == '__main__':
    description = 'Aligns, trims, and builds single gene trees from unaligned gene files.'
    usage = 'single_gene_tree_constructor.py -i path/to/input/ [OPTIONS]'
    parser, optional, required = help_formatter.initialize_argparse(name='single_gene_tree_constructor.py',
                                                                    desc=description,
                                                                    usage=usage)

    # Optional Arguments
    optional.add_argument('-t', '--threads', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Desired number of threads to be utilized.
                          Default: 1
                          """))
    optional.add_argument('--no_tree', action='store_true',
                          help=textwrap.dedent("""\
                              Do NOT build single gene trees.
                              Length filtration and trimmming only.
                              """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=True)

    # Parallelization of prepare_analyses function
    args.input = os.path.abspath(args.input)
    roots = [os.path.basename(file).split('.')[0] for file in os.listdir(args.input) if file.endswith('.fas')]
    args.output = os.path.abspath(args.output)
    mkdir_and_cd(args.output)
    mk_checkpoint_tmp(roots)
    checkpoints = get_checkpoints()
    file_count = len(roots)

    processes = args.threads
    if file_count < args.threads:
        processes = file_count

    prep_analyses = partial(prepare_analyses, checkpoints)
    with Pool(processes=processes) as p:
        all_checks = p.map(prep_analyses, roots)
