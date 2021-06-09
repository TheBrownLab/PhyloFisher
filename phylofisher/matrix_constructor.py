#!/usr/bin/env python
import csv
import os
import shutil
import subprocess
import sys
import textwrap
from collections import defaultdict
from glob import glob
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofisher import help_formatter


def bash(cmd):
    subprocess.run(cmd, executable='/bin/bash', shell=True)


def mk_dirs():
    """

    :return:
    """
    dirs = ['prequal', 'mafft', 'divvier', 'trimal']

    for my_dir in dirs:
        my_dir = f'{args.output}/{my_dir}'
        if os.path.isdir(my_dir) is False:
            os.mkdir(my_dir)


def delete_gaps_stars(gene, root):
    """Removes -'s and *'s from alignments"""
    file_name = f'{root}.aa'
    records = []
    for record in SeqIO.parse(gene, args.in_format):
        seq_string = str(record.seq).replace("-", "").replace("*", "")
        records.append(SeqRecord(Seq(seq_string),
                                 id='',
                                 name='',
                                 description=record.description))

    with open(file_name, 'w') as res:
        SeqIO.write(records, res, 'fasta')


def trim_and_align(gene):
    root = os.path.basename(gene).split('.')[0]
    # prequal
    os.chdir(f'{args.output}/prequal')
    delete_gaps_stars(gene, root)
    bash(f'prequal {root}.aa')
    os.chdir(args.output)

    # mafft
    os.chdir(f'{args.output}/mafft')
    bash(f'mafft --globalpair --maxiterate 1000 --unalignlevel 0.6 '
         f'--thread 1 {args.output}/prequal/{root}.aa.filtered > {root}.aln')
    os.chdir(args.output)

    # divvier
    os.chdir(f'{args.output}/divvier')
    bash(f'divvier -partial -mincol 4 -divvygap {args.output}/mafft/{root}.aln')
    shutil.move(f'{args.output}/mafft/{root}.aln.partial.fas', f'{args.output}/divvier/{root}.aln.partial.fas')
    shutil.move(f'{args.output}/mafft/{root}.aln.PP', f'{args.output}/divvier/{root}.aln.PP')
    os.chdir(args.output)

    # trimal
    os.chdir(f'{args.output}/trimal')
    bash(f'trimal -in {args.output}/divvier/{root}.aln.partial.fas -gt 0.80 -out {root}.gt80trimal.fas -fasta')
    os.chdir(args.output)


def parallelize(files):
    processes = args.threads
    if len(files) < args.threads:
        processes = len(files)

    with Pool(processes=processes) as p:
        all_checks = p.map(trim_and_align, files)


def parse_names(input_folder):
    """

    :param input_folder:
    :return:
    """
    name_set = set()
    if args.suffix:
        files = sorted(glob(f'{input_folder}/*{args.suffix}'))
    else:
        files = sorted(glob(f'{input_folder}/*'))
    for file in files:
        with open(file) as f:
            for record in SeqIO.parse(f, args.in_format):
                fname = record.description
                name = fname.split('_')[0]
                name_set.add(name)
    return files, sorted(list(name_set))


def stats(total_len, out_dict):
    """
    tools.py [OPTIONS] -i <input_dir> -m <metadata> {-n <gene_number> | -c <percent_complete>}

    :param total_len:
    :return:
    """
    with open(f'{args.output}/matrix_constructor_stats.tsv', 'w') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Taxon', 'PercentMissingData'])
        missing = []
        for record in SeqIO.parse(f'{args.output}/matrix.{out_dict[args.out_format.lower()]}', args.out_format.lower()):
            missing.append((record.name, (record.seq.count('-') / total_len) * 100))
        for org_missing in sorted(missing, key=lambda x: x[1], reverse=True):
            tsv_writer.writerow(list(org_missing))


if __name__ == '__main__':
    description = 'To trim align and concatenate orthologs into a super-matrix run:'
    parser, optional, required = help_formatter.initialize_argparse(name='matrix_constructor.py',
                                                                    desc=description,
                                                                    usage='matrix_constructor.py [OPTIONS] -i path/to/input/')

    # Optional Arguments
    optional.add_argument('-of', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Desired format of the output matrix.
                          Options: fasta, phylip (names truncated at 10 characters), 
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Format of the input files.
                          Options: fasta, phylip (names truncated at 10 characters), 
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))
    optional.add_argument('-c', '--concatenation_only', action='store_true', default=False,
                          help=textwrap.dedent("""\
                          Only concatenate alignments. Trimming and alignment are not performed automatically.
                          """))
    optional.add_argument('-t', '--threads', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Desired number of threads to be utilized.
                          Default: 1
                          """))

    # Changes help descriptions from the default input and output help descriptions
    in_help = 'Path to prep_final_dataset_<M.D.Y>'
    args = help_formatter.get_args(parser, optional, required, in_help=in_help)

    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)

    try:
        os.mkdir(args.output)
    except OSError:
        shutil.rmtree(args.output)
        os.mkdir(args.output)

    files, orgs = parse_names(args.input)
    if args.concatenation_only is False:
        os.chdir(args.output)
        mk_dirs()
        parallelize(files)
        files = sorted(glob(f'{args.output}/trimal/*'))

    total_len = 0
    res_dict = defaultdict(str)
    with open(f'{args.output}/indices.tsv', 'w') as outfile:
        outfile.write('Gene\tStart\tStop\n'
                      '')
        for file in files:
            gene = os.path.basename(file).split('.')[0]
            length = 0
            seq_dict = {}
            if args.concatenation_only:
                myformat = args.in_format
            else:
                myformat = 'fasta'
            for record in SeqIO.parse(file, myformat):
                length = len(record.seq)
                seq_dict[record.id.split('_')[0]] = str(record.seq)
            start_len = total_len + 1
            total_len += length
            outfile.write(f'{gene}\t{start_len}\t{total_len}\n')
            for org in orgs:
                if org in seq_dict:
                    res_dict[org] += seq_dict[org]
                else:
                    res_dict[org] += ('-' * length)

    # Accepted out formats with respective suffix
    out_dict = {'fasta':          'fas',
                'phylip':         'phy',
                'phylip-relaxed': 'phy',
                'nexus':          'nex'}

    # Creates SeqRecord iterator
    records = []
    for org, seq in res_dict.items():
        records.append(SeqRecord(Seq(seq),
                                 id=org,
                                 name='',
                                 description=''))

    # Writes to output matrix in user specified output
    if args.out_format.lower() in out_dict:
        with open(f'{args.output}/matrix.{out_dict[args.out_format.lower()]}', "w") as handle:
            SeqIO.write(records, handle, args.out_format.lower())
    else:
        sys.exit('Invalid Output Format')

    stats(total_len, out_dict)
