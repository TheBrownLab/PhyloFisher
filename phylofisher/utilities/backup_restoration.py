#!/usr/bin/env python
import os
import shutil
import tarfile
import textwrap
from datetime import datetime

from phylofisher import help_formatter


def list_backups():
    """
    List available backups in the database/backups dir
    :return:
    """
    print('Available Backups to Restore From:\n')
    print(f'Backup Num:\tDate:\t\tTime')
    for i, date_time in enumerate(backups):
        date = date_time.split('_')[0]
        time = date_time.split('_')[1]
        print(f'{i + 1}\t\t{date}\t{time}')


def restore():
    """
    Restores database from user specified backup
    :return:
    """

    if os.path.isdir(f'{args.database}/orthologs'):
        shutil.rmtree(f'{args.database}/orthologs')
    shutil.copytree(f'{backup_dir}/{backups[args.restore - 1]}/orthologs', f'{args.database}/orthologs')

    if os.path.isdir(f'{args.database}/paralogs'):
        shutil.rmtree(f'{args.database}/paralogs')
    shutil.copytree(f'{backup_dir}/{backups[args.restore - 1]}/paralogs', f'{args.database}/paralogs')

    if os.path.isfile(f'{args.database}/metadata.tsv'):
        os.remove(f'{args.database}/metadata.tsv')
    shutil.copy(f'{backup_dir}/{backups[args.restore - 1]}/metadata.tsv', f'{args.database}/metadata.tsv')

    if os.path.isfile(f'{args.database}/tree_colors.tsv'):
        os.remove(f'{args.database}/tree_colors.tsv')
    shutil.copy(f'{backup_dir}/{backups[args.restore - 1]}/tree_colors.tsv', f'{args.database}/tree_colors.csv')

    if os.path.isdir(f'{args.database}/proteomes'):
        shutil.rmtree(f'{args.database}/proteomes')
    shutil.copytree(f'{backup_dir}/{backups[args.restore - 1]}/proteomes', f'{args.database}/proteomes')


if __name__ == '__main__':
    desc = 'Utility to restore the PhyloFisher Database from backups'
    parser, optional, required = help_formatter.initialize_argparse(name='backup_restoration.py',
                                                                    desc=desc,
                                                                    usage="backup_restoration.py -d <db_dir>")

    required.add_argument('-d', '--database', type=str, metavar='<db_dir>', required=True,
                          help=textwrap.dedent("""\
                          Path to database directory.
                          """))

    optional.add_argument('--list_backups', action='store_true',
                          help=textwrap.dedent("""\
                          List available backups to restore from.
                          """))
    optional.add_argument('--restore', type=int, metavar='N',
                          help=textwrap.dedent("""\
                          Backup number to restore from. Backup number is obtained through the --list_backups option.
                          """))

    args = help_formatter.get_args(parser, optional, required, out_dir=False, pre_suf=False, inp_dir=False)

    backup_dir = f'{os.path.abspath(args.database)}/backups'
    backups = os.listdir(backup_dir)
    backups.sort(key=lambda date: datetime.strptime(date, '%d-%b-%Y_%H:%M:%S'), reverse=True)

    if args.list_backups:
        list_backups()

    if args.restore:
        restore()
