#!/usr/bin/env python
import os
import sys
import textwrap
from glob import glob
from phylofisher import help_formatter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from collections import defaultdict
import csv


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
            for line in f:
                if line.startswith('>'):
                    fname = line.split()[0][1:]
                    name = fname.split('_')[0]
                    name_set.add(name)
    return files, sorted(list(name_set))


def stats(total_len):
    """

    :param total_len:
    :return:
    """
    with open(f'{args.output}_forge_stats.tsv', 'w') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['org', 'missing[%]'])
        missing = []
        for record in SeqIO.parse(args.output, 'fasta'):
            missing.append((record.name, (record.seq.count('-') / total_len) * 100))
        for org_missing in sorted(missing, key=lambda x: x[1], reverse=True):
            tsv_writer.writerow(list(org_missing))


def main():
    """

    :return: NONE
    """
    input_folder = os.path.basename((args.input.strip('/')))
    files, orgs = parse_names(input_folder)
    total_len = 0
    res_dict = defaultdict(str)
    with open(f'{args.output}_indices.tsv', 'w') as outfile:
        for file in files:
            gene = os.path.basename(file).split('.')[0]
            length = 0
            seq_dict = {}
            for record in SeqIO.parse(file, 'fasta'):
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
    out_dict = {'fasta': 'fas',
                'phylip': 'phy',
                'phylip-relaxed': 'phy',
                'nexus': 'nex'}

    # Creates SeqRecord iterator
    records = []
    for org, seq in res_dict.items():
        records.append(SeqRecord(Seq(seq, IUPAC.protein),
                                 id=org,
                                 name='',
                                 description=''))

    # Writes to output matrix in user specified output
    if args.out_format.lower() in out_dict:
        with open(f'{args.output}', "w") as handle:
            SeqIO.write(records, handle, args.out_format.lower())
    else:
        sys.exit('Invalid Output Format')

    stats(total_len)


if __name__ == '__main__':
    parser, optional, required = help_formatter.initialize_argparse(name='forge.py',
                                                                    desc='Concatenates alignments into one'
                                                                         ' super-matrix.',
                                                                    usage='forge.py [OPTIONS] -i path/to/input/',
                                                                    dataset=False,
                                                                    input_meta=False)

    # Required Arguments
    optional.add_argument('-f', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Desired format of the output matrix.
                          Options: fasta, phylip (names truncated at 10 characters), 
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))

    # Changes help descriptions from the default input and output help descriptions
    in_help = 'Path to input directory containing alignments in FASTA format'
    args = help_formatter.get_args(parser, optional, required, in_help=in_help)

    main()
