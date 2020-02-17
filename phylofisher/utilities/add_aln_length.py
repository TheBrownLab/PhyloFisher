#!/usr/bin/env python
import sys
import os

def add_length(root):
    length = open(f'{root}.length').readline()
    os.rename(f'RAxML_bipartitions.{root}.tre', f'RAxML_bipartitions.{root}_{length}.tre')

add_length(sys.argv[1])