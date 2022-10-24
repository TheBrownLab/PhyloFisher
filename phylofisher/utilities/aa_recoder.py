#!/usr/bin/env python
from symbol import continue_stmt
import textwrap

from Bio import SeqIO
from Bio.Seq import Seq

from phylofisher import help_formatter

if __name__ == '__main__':
    description = 'Recodes input matrix based on SR4 amino acid classification.'
    parser, optional, required = help_formatter.initialize_argparse(name='aa_recoder.py',
                                                                    desc=description,
                                                                    usage='aa_recoder.py [OPTIONS] '
                                                                          '-i <matrix>'
                                                                          '-re <strategy>')
    # Required Arguments
    required.add_argument('-i', '--input', required=True, type=str, metavar='matrix',
                          help=textwrap.dedent("""\
                                    Path to input matrix for recoding.
                                    """))
                                    
    required.add_argument('-re', '--recoding_strat', metavar='strategy', type=str,
                          help=textwrap.dedent("""\
                                    Desired recoding strategy.
                                    Options: SR4, SR6, KGB6, D6, D9, D12, D15, or D18.
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

    # dictionary of recoding strategies
    key = dict(
    # SR4 class coding
        SR4 = {'A': ['A', 'G', 'N', 'P', 'S', 'T'],
            'C': ['C', 'H', 'W', 'Y'],
            'G': ['D', 'E', 'K', 'Q', 'R'],
            'T': ['F', 'I', 'L', 'M', 'V']
            },
            
        # SR6 class coding
        SR6 = {'0': ['A','P', 'S', 'T'],
            '1': ['D', 'E', 'N', 'G'],
            '2': ['Q', 'K', 'R'],
            '3': ['M', 'I', 'L', 'V'],
            '4': ['W', 'C'],
            '5': ['F', 'Y', 'H']
            },
            
        # KGB6 class coding
        KGB6 = {'0': ['A','P', 'G', 'S'],
            '1': ['D', 'E', 'N', 'Q', 'H', 'K', 'R', 'T'],
            '2': ['M', 'I', 'L'],
            '3': ['W'],
            '4': ['F', 'Y'],
            '5': ['C', 'V']
            },

        # Dayhoff 6-state coding
        D6 = {'0': ['A', 'G', 'P', 'S', 'T'],
            '1': ['D', 'E', 'N', 'Q'],
            '2': ['H', 'K', 'R'],
            '3': ['M', 'I', 'L', 'V'],
            '4': ['W', 'F', 'Y'],
            '5': ['C']
            },


        # Dayhoff 9-state coding
        D9 = {'0': ['D', 'E', 'H', 'N', 'Q'],
                '1': ['I', 'L','M','V'],
                '2': ['F', 'Y'],
                '3': ['A', 'S', 'T'],
                '4': ['K', 'R'],
                '5': ['G'],
                '6': ['P'],
                '7': ['C'],
                '8': ['W']
                },

        # Dayhoff 12-state coding
        D12 = {'0': ['D', 'E', 'Q'],
                '1': ['M', 'L', 'I', 'V'],
                '2': ['F', 'Y'],
                '3': ['K', 'H', 'R'],
                '4': ['G'],
                '5': ['A'],
                '6': ['P'],
                '7': ['S'],
                '8': ['T'],
                '9': ['N'],
                'A': ['W'],
                'B': ['C']
        },

        # Dayhoff 15-state coding
        D15 = {'0': ['D', 'E', 'Q'],
                '1': ['M', 'L'],
                '2': ['I', 'V'],
                '3': ['F', 'Y'],
                '4': ['G'],
                '5': ['A'],
                '6': ['P'],
                '7': ['S'],
                '8': ['T'],
                '9': ['N'],
                'A': ['K'],
                'B': ['H'],
                'C': ['R'],
                'D': ['W'],
                'E': ['C']
        },

        # Dayhoff 18-state coding
        D18 = {'0': ['M', 'L'],
                '1': ['F', 'Y'],
                '2': ['I'],
                '3': ['V'],
                '4': ['G'],
                '5': ['A'],
                '6': ['P'],
                '7': ['S'],
                '8': ['T'],
                '9': ['D'],
                'A': ['E'],
                'B': ['Q'],
                'C': ['N'],
                'D': ['K'],
                'E': ['H'],
                'F': ['R'],
                'G': ['W'],
                'H': ['C']
        }
    )

    out_dict = {'fasta'         : 'fas',
                'phylip'        : 'phy',
                'phylip-relaxed': 'phy',
                'nexus'         : 'nex'}

   
    #check input recoding strategy is valid
    if args.recoding_strat in key.keys():
        key = key[args.recoding_strat]
    else:
        error_msg = "Error: Invalid input recoding stragety \n Valid strateies include SR4, SR6, KGB6, Dayhoff6, Dayhoff9, Dayhoff12, Dayhoff15, and Dayhoff18"
        raise Exception(error_msg)
    
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
