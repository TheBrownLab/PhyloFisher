#!/usr/bin/env python
import os
import shutil
import textwrap
from statistics import mean

from Bio import SeqIO
from ete3 import Tree

from phylofisher import help_formatter


class Leaves:

    def __init__(self, tree):
        self.tree = tree
        self.leaf_names = self.tree.get_leaf_names()
        self.leaves = self.tree.get_leaves()

    def get_locations(self):
        locs = {}
        for name in self.leaf_names:
            locs[name] = self.tree & name
        return locs

    def get_mean_distance(self, leaf):
        locations = self.get_locations()
        distances = []
        for name in self.leaf_names:
            if name != leaf.name:
                distances.append(leaf.get_distance(locations[name]))
        sorted_dist = sorted(distances, reverse=True)
        return mean(sorted_dist[:   10])

    def org_speed(self):
        mean_distances = []
        for leaf in self.leaves:
            mean_distances.append((self.get_mean_distance(leaf), leaf.name))
        sorted_dist = sorted(mean_distances, reverse=True)
        return [i[1] for i in sorted_dist]


class Matrix:

    def __init__(self, matrix, format, ranked_orgs):
        self.matrix = matrix
        self.format = format
        self.ranked_orgs = ranked_orgs

    def fast_evol_taxa(self, step_size):
        steps = []
        for i in range(0, len(self.ranked_orgs), step_size):
            steps.append(self.ranked_orgs[0:i + step_size])
        return steps

    def generate_subset(self, output_folder, iterations, step_size):
        steps = self.fast_evol_taxa(step_size)
        os.mkdir(output_folder)
        for i in range(0, iterations):
            exclude = steps[i]

            # Original Code. Does NOT re-align and trim
            #
            # output_name = os.path.join(output_folder, f'step{i}')
            # with open(output_name, 'w') as res:
            #     for record in SeqIO.parse(self.matrix, self.format):
            #         if record.name not in exclude:
            #             res.write(f'>{record.name}\n{record.seq}\n')

            seq_files = [file for file in os.listdir(args.ortholog_files) if file.endswith('.fas')]

            os.mkdir(f'{output_folder}/tmp')

            for file in seq_files:
                with open(f'{args.ortholog_files}/{file}', 'r') as infile, open(f'{output_folder}/tmp/{file}', 'w') as outfile:
                    records = []
                    for record in SeqIO.parse(infile, 'fasta'):
                        if record.name not in exclude:
                            records.append(record)

                    SeqIO.write(records, outfile, 'fasta')

            os.system(f'matrix_constructor.py '
                      f'-i {output_folder}/tmp '
                      f'-o {output_folder}/step_{i} '
                      f'-of {args.out_format.lower()} '
                      f'-t {args.threads}')

            shutil.rmtree(f'{output_folder}/tmp')


if __name__ == '__main__':
    description = 'Removes the fastest evolving taxa, based on branch length.'
    parser, optional, required = help_formatter.initialize_argparse(name='fast_taxa_remover.py',
                                                                    desc=description,
                                                                    usage='fast_taxa_remover.py [OPTIONS] '
                                                                          '-tr <tree> '
                                                                          '-m <matrix> '
                                                                          '-i <num_of_iterations>')
    # Required Arguments
    required.add_argument('-m', '--matrix', required=True, type=str, metavar='matrix',
                          help=textwrap.dedent("""\
                          Path to input matrix.
                          """))
    required.add_argument('-tr', '--tree', required=True, type=str, metavar='tree',
                          help=textwrap.dedent("""\
                          Path to input tree.
                          """))
    required.add_argument('-i', '--iterations', metavar='N', type=int, required=True,
                          help=textwrap.dedent("""\
                          Number of iterations
                          """))
    required.add_argument('-or', '--ortholog_files', metavar='N', type=str, required=True,
                          help=textwrap.dedent("""\
                          Path to directory containing single gene files. This will be th path to 
                          prep_final_dataset_<M.D.Y> if used within the main PhyloFisher workflow.
                          """))

    # Optional Arguments
    optional.add_argument('-in_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Input matrix format if not FASTA.
                          Options: fasta, phylip (names truncated at 10 characters), 
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))
    optional.add_argument('-out_format', metavar='<format>', type=str, default='fasta',
                          help=textwrap.dedent("""\
                          Desired output format.
                          Options: fasta, phylip (names truncated at 10 characters),
                          phylip-relaxed (names are not truncated), or nexus.
                          Default: fasta
                          """))
    optional.add_argument('-s', '--step_size', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Number taxa removed per iteration.
                          Default: 1
                          """))
    optional.add_argument('-t', '--threads', metavar='N', type=int, default=1,
                          help=textwrap.dedent("""\
                          Desired number of threads to be utilized.
                          Default: 1
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

    tree = Tree(args.tree)
    x = Leaves(tree)
    m = Matrix(args.matrix, args.in_format, x.org_speed())
    m.generate_subset(args.output, args.iterations, args.step_size)
