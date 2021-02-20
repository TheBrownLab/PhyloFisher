#!/usr/bin/env python
import sys
import csv
import configparser
import textwrap
from pathlib import Path
from phylofisher import help_formatter


def check_input_meta(dataset_fold, input_meta):
    """
    Checks that all ids are unique in input_metadata and also between
    input_metadata and metadata. 
    Checks for illegal characters in input ids.
    """

    database_metadata_ids = set()
    with open(str(Path(dataset_fold, 'metadata.tsv')), 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            id_ = row[0]
            if id_ != 'Unique ID':
                database_metadata_ids.add(id_)

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
    optional.add_argument('--tree_colors', metavar='<omcl_data>',
                          help=textwrap.dedent("""\
                              Path to alternative single gene tree color configuration file.
                              """))

    config = configparser.ConfigParser()
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)


    check_input_meta(args.database_folder, args.input_file)

    if not args.orthomcl:
        args.orthomcl = str(Path(args.database_folder, 'orthomcl'))

    if not args.tree_colors:
        args.tree_colors = str(Path(args.database_folder, 'tree_colors.tsv'))

    with open('config.ini', 'w') as configfile:
        config['PATHS'] = {'database_folder': args.database_folder,
                           'input_file':     args.input_file,
                           'orthomcl':       args.orthomcl,
                           'color_conf':     args.tree_colors}
        config.write(configfile)
