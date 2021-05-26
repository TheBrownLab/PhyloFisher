#!/usr/bin/env python
import textwrap

from Bio import SeqIO
from Bio.Seq import Seq

from phylofisher import help_formatter

if __name__ == '__main__':
    description = 'Recodes input matrix based on SR4 amino acid classification.'
    parser, optional, required = help_formatter.initialize_argparse(name='SR4_class_recoder.py',
                                                                    desc=description,
                                                                    usage='SR4_class_recoder.py [OPTIONS] '
                                                                          '-i <matrix>')
    # Required Arguments
    required.add_argument('-i', '--input', required=True, type=str, metavar='matrix',
                          help=textwrap.dedent("""\
                                      Path to input matrix for recoding.
                                      """))
    # Optional Arguments
    optional.add_argument('-if', '--in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                                      Input file format if not FASTA.
                                      Options: fasta, phylip (names truncated at 10 characters), 
                                      phylip-relaxed (names are not truncated), or nexus.
                                      Default: fasta
                                      """))
    optional.add_argument('-of', '--out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                                  Desired output format.
                                  Options: fasta, phylip (names truncated at 10 characters), 
                                  phylip-relaxed (names are not truncated), or nexus.
                                  Default: fasta
                                  """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

    # SR4 class coding
    key = {'A': ['A', 'G', 'N', 'P', 'S', 'T'],
           'C': ['C', 'H', 'W', 'Y'],
           'G': ['D', 'E', 'K', 'Q', 'R'],
           'T': ['F', 'I', 'L', 'M', 'V']
           }

    out_dict = {'fasta'         : 'fas',
                'phylip'        : 'phy',
                'phylip-relaxed': 'phy',
                'nexus'         : 'nex'}

    # Opens input and output files
    with open(args.input, 'r') as infile, open(f'{args.output}.{out_dict[args.out_format]}', 'w') as outfile:
        all_records = []
        # SeqIO parses input file
        for record in SeqIO.parse(infile, format=args.in_format):
            # Iterates through the key dictionary's keys
            for nuc in key.keys():
                # Iterates through key's values
                for aa in key[nuc]:
                    # Replaces aa with nuc and X's with -'s
                    sequence = str(record.seq).replace(aa, nuc).replace("X", "-")
                    # Creates new SeqIO record
                    record.seq = Seq(sequence)
                    # Appends newly created SeqIO record to a list of SeqIO records to be used by SeqIO.write()
                    all_records.append(record)

        # Writes sequences to the outfile in user specified output format
        SeqIO.write(sequences=all_records, handle=outfile, format=args.out_format)
