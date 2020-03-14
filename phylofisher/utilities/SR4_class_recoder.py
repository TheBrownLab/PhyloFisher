import argparse
import textwrap
import phylofisher.help_formatter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein


def get_args():
    global args
    formatter = lambda prog: phylofisher.help_formatter.MyHelpFormatter(prog, max_help_position=100)
    parser = argparse.ArgumentParser(prog='SR4_class_recoder.py',
                                     description='some description',
                                     usage='some usage',
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent("""\
                                             additional information:
                                             """))
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    # Required Arguments
    required.add_argument('-i', '--input', required=True, type=str, metavar='',
                          help=textwrap.dedent("""\
                          Name of input alignment file.
                          """))
    required.add_argument('-o', '--output', type=str, required=True, metavar='',
                          help=textwrap.dedent("""\
                          Desired name of resulting output file.
                          """))
    # Optional Aruments
    optional.add_argument('-in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                                  Input file format if not FASTA.
                                  Options: fasta, phylip (names truncated at 10 characters), 
                                  phylip-relaxed (names are not truncated), or nexus.
                                  Default: fasta
                                  """))
    optional.add_argument('-out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                              Desired output format.
                              Options: fasta, phylip (names truncated at 10 characters), 
                              phylip-relaxed (names are not truncated), or nexus.
                              Default: fasta
                              """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                          Show this help message and exit.
                          """))
    parser._action_groups.append(optional)
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    # SR4 class coding
    key = {'A': ['A', 'G', 'N', 'P', 'S', 'T'],
           'C': ['C', 'H', 'W', 'Y'],
           'G': ['D', 'E', 'K', 'Q', 'R'],
           'T': ['F', 'I', 'L', 'M', 'V']
           }

    # Opens input and output files
    with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
        all_records = []
        # SeqIO parses input file
        for record in SeqIO.parse(infile, format=args.in_format, alphabet=generic_protein):
            # Iterates through the key dictionary's keys
            for nuc in key.keys():
                # Iterates through key's values
                for aa in key[nuc]:
                    # Replaces aa with nuc and X's with -'s
                    sequence = str(record.seq).replace(aa, nuc).replace("X", "-")
                    # Creates new SeqIO record
                    record.seq = Seq(sequence, alphabet=generic_protein)
                    # Appends newly created SeqIO record to a list of SeqIO records to be used by SeqIO.write()
                    all_records.append(record)

        # Writes sequences to the outfile in user specified output format
        SeqIO.write(sequences=all_records, handle=outfile, format=args.out_format)
