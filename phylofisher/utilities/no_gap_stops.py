#!/usr/bin/env python
import sys
from Bio import SeqIO

def delete_gaps_stars(file):
    file_name = f'{file.split(".")[0]}.aa'
    with open(file_name, 'w') as res:
        for record in SeqIO.parse(file, 'fasta'):
            res.write(f'>{record.name}\n{str(record.seq).replace("-", "").replace("*", "")}\n')

delete_gaps_stars(sys.argv[1])