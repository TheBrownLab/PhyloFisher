#!/usr/bin/env python
import textwrap
import os
from ete3 import Tree
from phylofisher import help_formatter


def get_names(metadata):
    '''
    Parses metadata file to get short and long names for taxa

    :param metadata: path to metadata file
    :type metadata: str
    :return: dictionary of short and long names
    :rtype: dict
    '''
    name_dict = dict()
    with open(metadata, 'r') as infile:
        infile.readline()
        for line in infile:
            unique_id, long_name = line.split('\t')[0], line.split('\t')[1]
            if long_name in name_dict.values():
                long_name = long_name + '_' + unique_id
            
            name_dict[unique_id] = long_name
    
    return name_dict


if __name__ == '__main__':
    description = 'Renames leaves in a tree to match the names in the database.'
    parser, optional, required = help_formatter.initialize_argparse(name='leaf_renamer.py',
                                                                    desc=description,
                                                                    usage='leaf_renamer.py '
                                                                          '-tr <tree> '
                                                                          '-d <database> ')
   
    # Required Arguments
    required.add_argument('-d', '--database', required=True, type=str, metavar='database',
                          help=textwrap.dedent("""\
                          Path to PhyloFisher database.
                          """))
    required.add_argument('-tr', '--tree', required=True, type=str, metavar='tree',
                          help=textwrap.dedent("""\
                          Path to input tree.
                          """))
    required.add_argument('-o', '--output', required=True, type=str, metavar='tree',
                          help=textwrap.dedent("""\
                          Output tree name.
                          """))
    
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False, out_dir=False)

    metadata = str(os.path.join(args.database, 'metadata.tsv'))
    tree = Tree(args.tree, format=1)
    name_dict = get_names(metadata)
    
    # Update leaf names
    for node in tree.traverse():
        if node.is_leaf():
            node.name = name_dict[node.name]

    tree.write(outfile=args.output, format=1)