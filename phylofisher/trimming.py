#!/usr/bin/env python
import os
import argparse
from glob import glob
from pathlib import Path



def prepare_analyses(dataset, threads):
    root = dataset.split('/')[-1].split('.')[0]
    command = ""

    command += f'no_gap_stops.py {dataset} && '

    command += f'mafft --globalpair --maxiterate 1000 --unalignlevel 0.6' \
        f' --threads {threads} {dataset}.aa > {root}.aln && '

    divvier = str(Path(dfo, 'lib/Divvier/divvier'))
    command += f'{divvier} -mincol 4 -divvygap {root}.aln && '

    command += f'trimal -in {root}.aln.divvy.fas -gt 0.2 -out {root}.trimal\n'

    return command


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for deleting orgs/taxonomic'
                                                 'groups from the dataset', usage="purge.py [OPTIONS]")
    parser.add_argument('-d', '--dataset_folder')
    parser.add_argument('-i', '--input_folder', help='Short names of organisms for deletion: Org1,Org2,Org3')
    parser.add_argument('-s', '--suffix', default='.fas')
    parser.add_argument('-c', '--use_config', action='store_true')
    args = parser.parse_args()

    if args.use_config:
        config = configparser.ConfigParser()
        config.read('config.ini')
        dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
    else:
        dfo = args.dataset_folder

    os.chdir(args.input_folder)
    with open('commands.txt', 'w') as res:
        for file in glob(f'*{args.suffix}'):
            line = prepare_analyses(file, 5)
            res.write(line)
