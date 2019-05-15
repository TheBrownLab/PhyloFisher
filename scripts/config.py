#!/usr/bin/env python
import os
import configparser
import argparse
from pathlib import Path

config = configparser.ConfigParser()
parser = argparse.ArgumentParser(description='Script for parsing fun.',
                                 usage="config.py -d <dataset_folder> -i <input_file> [OPTIONS]")
parser.add_argument('-d', '--dataset_folder', required=True)
parser.add_argument('-i', '--input_file', required=True)
parser.add_argument('-o', '--orthomcl', help='path to orthomcl if not in dataset_folder')
# parser.add_argument('-s', '--stable_dataset', help='Dataset for long-term usage.')
args = parser.parse_args()

if not args.orthomcl:
    args.orthomcl = str(Path(args.dataset_folder, 'orthomcl'))

with open('config.ini', 'w') as configfile:
    config['PATHS'] = {'dataset_folder': args.dataset_folder,
                       'input_file': args.input_file,
                       'orthomcl': args.orthomcl}
    config.write(configfile)
