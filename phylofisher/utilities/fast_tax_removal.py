#!/usr/bin/env python
import os
from ete3 import Tree
from statistics import mean
from Bio import SeqIO

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
        print(chunks)
        return chunks


    def generate_subset(self, output_name, iterations, chunk_size=1):
        chunks = self.fast_evol_taxa(chunk_size)
        for i in range(0, iterations):
            exclude = chunks[i]
            with open(f'{output_name}_chunk{i}', 'w') as res:
                for record in SeqIO.parse(self.matrix, self.format):
                    if record.name not in exclude:
                        res.write(f'>{record.name}\n{record.seq}\n')



data_folder = '/home/david/MsPhylo/test/new/ForDavid_new/ForDavid/FastSiteRemoval'
tree_file = os.path.join(data_folder, 'Bordor.351.64.3-7-2016.dat.treefile')
matrix = os.path.join(data_folder, 'Bordor.351.64.3-7-2016.dat')

tree = Tree(tree_file)
x = Leaves(tree)
m = Matrix(matrix, 'phylip', x.org_speed())
m.generate_subset('jahody', 3)