from glob import glob
import argparse
from ete3 import Tree
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
# plt.style.use('ggplot')
plt.rcParams["figure.figsize"] = (10,5)

def bipartitions(tree):
    bipar = set()
    all_ = frozenset(tree.get_leaf_names())
    for node in tree.traverse('preorder'):
        if node.is_leaf() is False:
            taxa = frozenset(node.get_leaf_names())
            if len(taxa) < (len(all_) - len(taxa)):
                if len(taxa) > 1:
                    bipar.add(taxa)
            else:
                rest = taxa.symmetric_difference(all_)
                if len(rest) > 1:
                    bipar.add(rest)
    return bipar

def support(trees):
    bootstrap = []
    n_trees = 0
    for line in open(trees):
        tree = Tree(line)
        n_trees += 1
        bootstrap = bootstrap + list(bipartitions(tree))
    norm_bootstrap = {}
    counted = dict(Counter(bootstrap))
    for key, value in counted.items():
        norm_bootstrap[key] = value/n_trees
    return norm_bootstrap

def get_support(group, supp_dict):
    query = frozenset(group)
    if query in supp_dict:
        return supp_dict[query]
    else:
        return 0

def parse_groups(input_file):
    query_dict = {}
    for line in open(input_file):
        group, orgs = line.split(':')
        query_dict[group] = [org.strip() for org in orgs.split(',')]
    return query_dict.items()

def file_to_series(file):
    sup_dict = support(file)
    group_sup = {}
    for group, orgs in queries:
        group_sup[group] = get_support(orgs, sup_dict)
    s = pd.Series(group_sup)
    return s

def parse_bss():
    columns = []
    n = 0
    for line in open(args.bs_files):
        column = file_to_series(line.strip())
        if n == 0:
            column.name = 'Full dataset'
        else:
            column.name = f'Step {n}'
        columns.append(column)
        n += 1
    return columns

def main():
    columns = parse_bss()
    df = pd.DataFrame(columns)
    # df.sort_index(inplace=True)
    df.plot()
    plt.legend(loc=0, prop={'size': 6})
    plt.xticks(range(len(df.index)), list(df.index))
    plt.tight_layout()
    plt.savefig("test_out.pdf")
    df.to_csv("test_out.csv")
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='some description', usage="blabla")
    parser.add_argument('-b', '--bs_files')
    parser.add_argument('-g', '--groups')
    args = parser.parse_args()
    queries = parse_groups(args.groups)
    main()