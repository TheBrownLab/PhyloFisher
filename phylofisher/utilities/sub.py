#!/usr/bin/env python
import argparse

def make_sh(line, header):
    with open(f'{line.split()[1].split(".")[0]}.sh', 'w') as res:
        res.write(f'{header}\n{line}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for deleting orgs/taxonomic'
                                                 'groups from the dataset', usage="purge.py [OPTIONS]")
    parser.add_argument('-e', '--header')
    parser.add_argument('-c', '--commands')
    args = parser.parse_args()

    header = open(args.header).read()
    for line in open(args.commands):
        make_sh(line, header)