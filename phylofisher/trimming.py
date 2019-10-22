#!/usr/bin/env python
import os
import argparse
from glob import glob
from pathlib import Path
import configparser

#TODO fucking cleaning

def prepare_analyses(dataset, threads):
    root = dataset.split('/')[-1].split('.')[0]
    command = ""

    command += f'no_gap_stops.py {dataset} && '

    preq = str(Path(dfo, 'lib/prequal/prequal'))

    command += f'{preq} {root}.aa && '

    command += f' mafft --globalpair --maxiterate 1000 --unalignlevel 0.6' \
        f' --thread {threads} {root}.aa.filtered > {root}.aln && '

    divvier = str(Path(dfo, 'lib/Divvier/divvier'))
    command += f'{divvier} -mincol 4 -divvygap {root}.pre_aln && '

    command += f'pre_trimal.py {root}.pre_aln'

    command += f'trimal -in {root}.aln.divvy.fas -gt 0.2 -out {root}.trimal && '

    command += f'len_filter2.py -i {root}.trimal -t 0.5 && '

    command += f'mafft --globalpair --maxiterate 1000 --unalignlevel 0.6' \
        f' --thread {threads} {root}.len > {root}.aln2 && '

    command += f'{divvier} -mincol 4 -divvygap {root}.aln2 && '

    command += f'trimal -in {root}.aln2.divvy.fas -gt 0.01 -out {root}.final && '

    command += f'raxmlHPC-PTHREADS-AVX2 -T {threads} -m PROTGAMMALG4XF -f a -s {root}.final' \
        f' -n {root}.tre -x 123 -N 100 -p 12345 && '

    command += f'add_aln_length.py {root}\n'

    return command


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for deleting orgs/taxonomic'
                                                 'groups from the dataset', usage="purge.py [OPTIONS]")
    parser.add_argument('-d', '--dataset_folder')
    parser.add_argument('-i', '--input_folder', help='Short names of organisms for deletion: Org1,Org2,Org3')
    parser.add_argument('-s', '--suffix', default='.fas')
    parser.add_argument('-c', '--use_config', action='store_true')
    parser.add_argument('-t', '--threads', type=int, default=1)
    args = parser.parse_args()

    if args.use_config:
        config = configparser.ConfigParser()
        config.read('config.ini')
        dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
    else:
        dfo = str(Path(args.dataset_folder).resolve())

    os.chdir(args.input_folder)
    with open('commands.txt', 'w') as res:
        for file in glob(f'*{args.suffix}'):
            line = prepare_analyses(file, args.threads)
            res.write(line)
