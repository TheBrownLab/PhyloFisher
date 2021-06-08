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
    print(f'Backup Num:\tDate:\t\tTime:')
    for i, date_time in enumerate(backups):
        date_time_str = datetime.strftime(date_time, "%d-%b-%Y_%H-%M-%S")
        date = date_time_str.split('_')[0]
        time = date_time_str.split('_')[1]
        print(f'{i + 1}\t\t{date.replace("-"," ")}\t{time.replace("-",":")}')


def restore():
    """
    Restores database from user specified backup
    :return:
    """
    backup_dir = f'{args.database}/backups/{datetime.strftime(backups[args.restore - 1], "%d-%b-%Y_%H-%M-%S")}'

    if os.path.isdir(f'{args.database}/orthologs'):
        shutil.rmtree(f'{args.database}/orthologs')
    shutil.copytree(f'{backup_dir}/orthologs', f'{args.database}/orthologs')

    if os.path.isdir(f'{args.database}/paralogs'):
        shutil.rmtree(f'{args.database}/paralogs')
    shutil.copytree(f'{backup_dir}/paralogs', f'{args.database}/paralogs')

    if os.path.isfile(f'{args.database}/metadata.tsv'):
        os.remove(f'{args.database}/metadata.tsv')
    shutil.copy(f'{backup_dir}/metadata.tsv', f'{args.database}/metadata.tsv')

    if os.path.isfile(f'{args.database}/tree_colors.tsv'):
        os.remove(f'{args.database}/tree_colors.tsv')
    shutil.copy(f'{backup_dir}/tree_colors.tsv', f'{args.database}/tree_colors.tsv')

    if os.path.isdir(f'{args.database}/proteomes'):
        shutil.rmtree(f'{args.database}/proteomes')
    shutil.copytree(f'{backup_dir}/proteomes', f'{args.database}/proteomes')


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
    backups = [datetime.strptime(date, '%d-%b-%Y_%H-%M-%S') for date in backups]
    backups.sort()

    if args.list_backups:
        list_backups()

    if args.restore:
        restore()
