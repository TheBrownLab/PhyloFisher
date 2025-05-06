#!/usr/bin/env python
import sys
import csv
import configparser
import textwrap
import sqlite3
from pathlib import Path
from phylofisher import help_formatter


def check_input_meta(input_meta):
    '''
    Checks the input metadata file for unique IDs and illegal characters.

    :param input_meta: path to the input metadata file in tsv format
    :type input_meta: str
    '''

    cursor.execute('SELECT short_name FROM metadata')
    database_metadata_ids = set([x[0] for x in cursor.fetchall()])

    input_ids = set()
    with open(input_meta, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            id_ = row[2]
            if id_ == 'Unique ID':
                continue
            
            #check that ID is unique
            if id_ in input_ids:
                print("ERROR: ", id_, "is not a unique ID")
                sys.exit()

            if id_ in database_metadata_ids:
                print("ERROR: ", id_, "already in the database")
                sys.exit()

            #check illegal characters
            for ch in ["_", '@', '..', '*', ' ']:
                if ch in id_:
                    print(f"ERROR: illegal character {ch} in {id_}")
                    sys.exit()

            input_ids.add(id_)


if __name__ == '__main__':
    description = 'Script for the analysis folder configuration.'
    usage = 'config.py -d <database_folder> -i <input_file.tsv> [OPTIONS]'
    parser, optional, required = help_formatter.initialize_argparse(name='config.py',
                                                                    desc=description,
                                                                    usage=usage)

    # Required Arguments
    required.add_argument('-d', '--database_folder', required=True, metavar='<database_dir>',
                          help=textwrap.dedent("""\
                          Path to directory containing the PhyloFisher database.
                          """))
    required.add_argument('-i', '--input_file', required=True, metavar='<input.tsv>',
                          help=textwrap.dedent("""\
                          Path to input metadata in tsv format.
                          """))
    # Optional Arguments
    optional.add_argument('--orthomcl', metavar='<omcl_data>',
                          help=textwrap.dedent("""\
                          Path to orthomcl if NOT in PhyloFisherDatabase_v1.0/database/orthomcl
                          """))

    config = configparser.ConfigParser()
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)


    conn = sqlite3.connect(str(Path(args.database_folder, 'phylofisher.db')))
    cursor = conn.cursor()
    check_input_meta(conn, args.input_file)

    if not args.orthomcl:
        args.orthomcl = str(Path(args.database_folder, 'orthomcl'))


    with open('config.ini', 'w') as configfile:
        config['PATHS'] = {'database_folder': args.database_folder,
                           'input_file':     args.input_file,
                           'orthomcl':       args.orthomcl}
        config.write(configfile)
