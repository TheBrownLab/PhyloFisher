#!/usr/bin/env python
import sys
from Bio import SeqIO

def x_to_dash(file):
    file_name = f'{file.split(".")[0]}.pre_trimal'
    with open(file_name, 'w') as res:
        for record in SeqIO.parse(file, 'fasta'):
            res.write(f'>{record.name}\n{str(record.seq).replace("X", "-")}\n')

x_to_dash(sys.argv[1])