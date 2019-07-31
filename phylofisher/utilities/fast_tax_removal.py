#!/usr/bin/env python
import os
from ete3 import Tree
from statistics import mean
from Bio import SeqIO
import argparse

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
        return mean(sorted_dist[:10])


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


    def fast_evol_taxa(self, chunk_size):
        chunks = []
        for i in range(0, len(self.ranked_orgs), chunk_size):
            chunks.append(self.ranked_orgs[0:i+chunk_size])
        return chunks


    def generate_subset(self, output_folder, iterations, chunk_size):
        chunks = self.fast_evol_taxa(chunk_size)
        os.mkdir(output_folder)
        for i in range(0, iterations):
            exclude = chunks[i]
            output_name = os.path.join(output_folder, f'chunk{i}')
            with open(output_name, 'w') as res:
                for record in SeqIO.parse(self.matrix, self.format):
                    if record.name not in exclude:
                        res.write(f'>{record.name}\n{record.seq}\n')


def main():
    parser = argparse.ArgumentParser(description='Fast taxon removal', usage="fast_tax_removal.py [OPTIONS]")
    parser.add_argument('-t', '--tree', required=True)
    parser.add_argument('-m', '--matrix', required=True)
    parser.add_argument('-f', '--format', default='fasta', help='format of your matrix [default: fasta]')
    parser.add_argument('-o', '--output_folder')
    parser.add_argument('-i', '--iterations', required=True, type=int)
    parser.add_argument('-c', '--chunk_size', default=1, type=int)
    args = parser.parse_args()

    tree = Tree(args.tree)
    x = Leaves(tree)
    m = Matrix(args.matrix, args.format, x.org_speed())
    m.generate_subset(args.output_folder, args.iterations, args.chunk_size)


if __name__ == '__main__':
    main()