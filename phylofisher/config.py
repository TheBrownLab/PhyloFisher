#!/usr/bin/env python
import configparser
import argparse
from pathlib import Path

config = configparser.ConfigParser()
parser = argparse.ArgumentParser(description='Script for parsing fun.',
                                 usage="config.py -d <dataset_folder> -i <input_file> [OPTIONS]")
parser.add_argument('-d', '--dataset_folder', required=True)
parser.add_argument('-i', '--input_file', required=True)
parser.add_argument('--bmge', help='path to BMGE jar file')
parser.add_argument('--orthomcl', help='path to orthomcl if not in dataset_folder')
args = parser.parse_args()

if not args.orthomcl:
    args.orthomcl = str(Path(args.dataset_folder, 'orthomcl'))

with open('config.ini', 'w') as configfile:
    config['PATHS'] = {'dataset_folder': args.dataset_folder,
                       'input_file': args.input_file,
                       'orthomcl': args.orthomcl,
                       'bmge': args.bmge}
    config.write(configfile)
