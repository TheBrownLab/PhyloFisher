#!/usr/bin/env python
import sys
import glob
import os
import argparse
from collections import defaultdict
from collections import Counter
from multiprocessing import Pool
from pathlib import Path
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, PieChartFace
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.style.use('ggplot')



def parse_metadata():
    metadata = {}
    tax_col = {}
    for line_ in open(str(Path(dfo, 'metadata.tsv'))):
        if 'Full Name' not in line_:
            sline = line_.split('\t')
            tax = sline[0].strip()
            group = sline[2].strip()
            col = sline[3].strip()
            full = sline[1].strip()
            metadata[tax] = {'group': group, 'col': col, 'full': full}
            tax_col[group] = col
    for line in open(multi_input):
        if "FILE_NAME" not in line:
            metadata_input = line.split('\t')
            tax = metadata_input[2].strip().split('_')[0]
            group = metadata_input[3].strip()
            full = metadata_input[5].strip()
            metadata[tax] = {'group': group, 'col': "white", 'full': full}
    print(metadata)
    return metadata, tax_col

def suspicious_clades(tree):
    t = Tree(tree)
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)

    supported_clades = []
    for node in t.traverse('preorder'):
        if (node.is_root() is False) and (node.is_leaf() is False):
            if node.support >= 70 and (len(node) < (len(t)-len(node))):
                clade = node.get_leaf_names()
                if len(clade) > 1:
                    supported_clades.append(clade)
    suspicious = []
    for clade in supported_clades:
        groups = set()
        for org in clade:
            org = org.split('_')[0]
            groups.add(metadata[org]['group'])
        if len(groups) > 1:
            suspicious.append(clade)
    return tree, suspicious


def get_best_candidates(tree_file):
    t = Tree(tree_file)
    top_rank = defaultdict(dict)
    for node in t.traverse('preorder'):
        if node.is_leaf():
            if node.name.count('_') == 4:
                org,__, _, rank, _ = node.name.split('_')
                rank = int(rank[1:])
                if org not in top_rank:
                    top_rank[org]['rank'] = rank
                    top_rank[org]['candidate'] = node.name
                else:
                    if rank < top_rank[org]['rank']:
                        top_rank[org]['rank'] = rank
                        top_rank[org]['candidate'] = node.name
    top_seqs = set()
    for org_ in top_rank.values():
        top_seqs.add(org_['candidate'])
    return top_seqs


def fucking_shit(names):
    for name in names:
        if "GEN" in name:
            if str(names).count(name.replace('GEN', '')) == 2:
                return True
    return False


def tree_to_pdf(tree_file):
    top_ranked = get_best_candidates(tree_file)
    output_base = f"{output_folder}/{tree_file.split('.')[1]}.pdf"
    table = open(f"{output_folder}/{tree_file.split('.')[1]}.tsv", 'w')
    t = Tree(tree_file)
    ts = TreeStyle()
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    sus_clades = 0

    to_delete = ["PhysTSA", "Acancast","Thalpseu","Tricvagi","Emilhuxl",
                   "Dictdisc","Plasbras","Sacccere","Ectosili",
                   "Chlarein","Planfung","Monoexil","Plasfalc",
                   "Bigenata","Paratetr","Porppurp","Spirsalm",
                   "Phytpara","Batrdend","Tetrther","Cracchor",
                   "Kipfbial","Naeggrub","Thectrah","Vitrbras",
                   "Charbrau","Entahist","Monobrev",
                   "Symbmicr","Goniavon","Homosapi","Nemavect",
                   "Ichtmult","Trypbruc","Guilthet","Capsowcz",
                   "Oxyrmari", 'Oxyttrif', 'CracchorGEN', 'Colpangu2', 'Pyroyezo']

    for node in t.traverse('preorder'):
        node_style = NodeStyle()
        node_style['vt_line_width'] = 3
        node_style['hz_line_width'] = 3
        node_style['vt_line_type'] = 0
        node_style['hz_line_type'] = 0
        if node.support >= 70:
            node_style['shape'] = 'circle'
            node_style['size'] = 12
            node_style['fgcolor'] = 'black'

        if node.is_root() is False:
            if node.is_leaf() is False:
                supp = TextFace(f'{int(node.support)}', fsize=8)
                if node.support >= 70:
                    supp.bold = True
                    taxons = set()
                    taxons_list = []
                    orgs = node.get_leaf_names()
                    if len(orgs) > 1:
                        for org in orgs:
                            org = org.split('_')[0]
                            taxons.add(metadata[org]['group'])
                            taxons_list.append(metadata[org]['group'])
                    if len(taxons) > 1 and (len(node) < (len(t)-len(node))):
                        node_style['shape'] = 'sphere'
                        node_style['fgcolor'] = 'red'
                        node_style['bgcolor'] = 'Silver'
                        sus_clades += 1
                else:
                    supp.fsize = 7
                node.add_face(supp, column=0, position="branch-bottom")
            else:
                empty_face = TextFace("\t"*20)
                node.add_face(empty_face, column=2, position = "aligned")
                node.add_face(empty_face, column=3, position="aligned")
                org = node.name
                if org.count('_') == 4:
                    quality = f'{org.split("_")[-3]}_{org.split("_")[-2]}_{org.split("_")[-1]}'
                    org = org.split('_')[0]
                    group = metadata[org]['group']
                    color = metadata[org]['col']
                    node_style["bgcolor"] = color
                    node.name = f'{node.name}'
                    if group in tax_col:
                        tax_name = TextFace(f'[{group}]', fgcolor=tax_col[group], bold=True)
                    else:
                        tax_name = TextFace(f'[{group}]', bold=True)
                    node.add_face(tax_name, column=1, position = "aligned") #TODO FUUUUCK ME
                    if org in to_delete:
                        table.write(f'{metadata[org]["full"]}_{quality}@{org}\t{group}\td*\n')
                        deletef = TextFace(f'{metadata[org]["full"]}_{quality}@{org}', fgcolor='red')
                        node.name = ''
                        node.add_face(deletef, column=0)
                    elif node.name in contaminations:
                        table.write(f'{metadata[org]["full"]}_{quality}@{org}\t{group}\td*\n')
                        deletef = TextFace(f'{metadata[org]["full"]}_{quality}@{org}', fgcolor='red')
                        node.name = ''
                        node.add_face(deletef, column=0)
                    elif node.name in top_ranked:
                        tname = TextFace(f'{metadata[org]["full"]}_{quality}@{org}', bold=True)
                        table.write(f'{metadata[org]["full"]}_{quality}@{org}\t{group}\to*\n')
                        node.name = ''
                        node.add_face(tname, column=0, position='branch-right')
                    else:
                        table.write(f'{metadata[org]["full"]}_{quality}@{org}\t{group}\tp*\n')
                        node.name = f'{metadata[org]["full"]}_{quality}@{org}'
                else:
                    org, length = org.split('_')
                    group = f"[{metadata[org]['group']}]"
                    gface = TextFace(f'{group}') #TODO do not touch me pleeeease
                    if org in to_delete:
                        gface = TextFace(f'{group}', fgcolor=tax_col[group[1:-1]], bold=True)
                        table.write(f'{metadata[org]["full"]}_{length}@{org}\t{group[1:-1]}\td*\n')
                        deletef = TextFace(f'{metadata[org]["full"]}_{length}@{org}', fgcolor='red')
                        node.name = ''
                        node.add_face(deletef, column=0)
                        node.add_face(gface, column=1, position="aligned")
                    else:
                        color = metadata[org]['col']
                        node_style["bgcolor"] = color
                        table.write(f'{metadata[org]["full"]}_{length}@{org}\t{group[1:-1]}\to\n')
                        node.name = f'{metadata[org]["full"]}_{length}@{org}'
                        node.add_face(gface, column=1, position="aligned")

            node.set_style(node_style)

    title_face = TextFace(f'<{output_base.split("/")[-1]}, {sus_clades} suspicious clades>',  bold=True)
    ts.title.add_face(title_face, column=1)
    t.render(output_base + '_tree.svg', tree_style=ts)
    table.close()


def trees_plus_table(trees):
    # with Pool(processes=threads) as pool:
    #     pool.map(tree_to_pdf, trees)
    with Pool(processes=threads) as pool:
        suspicious = list(pool.map(suspicious_clades, trees))
        return suspicious


def nonredundant(result_clades):
    nonredundant = []
    for clades in result_clades:
        if clades:
            sorted_clades = sorted(clades, key=lambda clade: len(clade))
            sorted_sets = [set(clade) for clade in sorted_clades]
            for i, set_ in enumerate(sorted_sets):
                redundant = False
                if i < len(sorted_sets):
                    for other in sorted_sets[i+1:]:
                        if set_.issubset(other):
                            redundant = True
                if redundant is False:
                    nonredundant.append(set_)
    return nonredundant


def problematic(nonredundant_clades_):
    problematic_orgs_ = defaultdict(list)
    contaminations = set()
    for clade in nonredundant_clades_:
        taxons = []
        for org in clade:
            org = org.split('_')[0]
            taxons.append(metadata[org]["group"])
        group_perc = {}
        for tax in set(taxons):
            group_perc[tax] = (taxons.count(tax)/len(taxons)) * 100
        max_perc = max(group_perc, key=group_perc.get)
        for seq in clade:
            org = seq.split('_')[0]

            if group_perc[metadata[org]["group"]] == 50:
                problematic_orgs_[org].append('unspecified')
            else:
                if org == 'AcanthFB' and max_perc == 'Discoba':
                    contaminations.add(seq)

                elif org == 'Aurisoli' and max_perc == 'Obazoa':
                    contaminations.add(seq)

                elif org == 'ChoanoFB' and max_perc in ['Obazoa', 'Discoba']:
                    contaminations.add(seq)

                elif org == 'Colpangu1' and max_perc in ['Stramenopiles', 'Discoba']:
                    contaminations.add(seq)

                elif org == 'Lapoguse' and max_perc in ['Discoba']:
                    contaminations.add(seq)

                elif org == 'Leptvora' and max_perc in ['Chloroplastida', 'Obazoa']:
                    contaminations.add(seq)

                elif org == 'Mincchit' and max_perc in ['Stramenopiles', 'Chloroplastida']:
                    contaminations.add(seq)

                elif org == 'Oxyrmar2' and max_perc in ['Stramenopiles']:
                    contaminations.add(seq)

                elif org == 'Platmacr' and max_perc == 'Discoba':
                    contaminations.add(seq)

                elif org == 'Pyroyezo' and max_perc in ['Alveolata']:
                        contaminations.add(seq)

                elif org == 'Sori0191' and max_perc in ['Alveolata', 'Obazoa']:
                        contaminations.add(seq)

                elif org == 'Stichzan' and max_perc in ['CRuMs']:
                        contaminations.add(seq)

                elif org == 'SymbmicrGEN' and max_perc in ['Obazoa', 'Rhizaria']:
                        contaminations.add(seq)

                elif org == 'Syssmult' and max_perc in ['Stramenopiles']:
                        contaminations.add(seq)

                elif org == 'TeloneP1' and max_perc in ['Discoba']:
                        contaminations.add(seq)

                elif org == 'TeloneP2' and max_perc in ['Discoba']:
                        contaminations.add(seq)

                elif org == 'Telosubt' and max_perc in ['Discoba']:
                        contaminations.add(seq)

                problematic_orgs_[org].append(max_perc)


    return problematic_orgs_, contaminations


def plot_problematic(problematic_orgs):
    max_ = 0
    print(problematic_orgs)
    for org, groups in sorted(problematic_orgs.items()):
        counted_groups = dict(Counter(groups))
        group_total = 0
        for group in counted_groups.values():
            group_total += group
        if group_total > max_:
            max_ = group_total
    print(max_)
    with PdfPages('problematic_orgs.pdf') as pdf:
        for org, groups in sorted(problematic_orgs.items()):
            counted_groups = dict(Counter(groups))
            colors = []
            for group in counted_groups.keys():
                if group == 'unspecified':
                    colors.append('black')
                else:
                    colors.append(tax_col.get(group, "white"))
            barplot = plt.bar(counted_groups.keys(), counted_groups.values(), color=colors)
            plt.tight_layout()
            plt.ylim(0, max_)
            plt.title(f'{org}({metadata[org]["group"]})')
            plt.xticks(rotation='vertical')
            plt.yticks(fontsize=4)
            pdf.savefig(bbox_inches='tight')
            plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='some description', usage="blabla")
    parser.add_argument('-tf', '--trees_folder')
    parser.add_argument('-o', '--output_folder')
    parser.add_argument('-df', '--dataset_folder')
    parser.add_argument('-i', '--infile')
    args = parser.parse_args()
    trees_folder = args.trees_folder
    output_folder = args.output_folder
    dfo = args.dataset_folder
    multi_input = os.path.abspath(args.infile)
    os.mkdir(output_folder)
    trees = glob.glob(f"{trees_folder}/*tre*")
    number_of_genes = len(trees)
    metadata, tax_col = parse_metadata()
    threads = 1
    suspicious = trees_plus_table(trees)
    suspicious = sorted(suspicious)

    with open('suspicious.txt', 'w') as res:
        for file, clades in suspicious:
            if clades:
                res.write(f'{file}\n========================\n')
                for clade_ in clades:
                    res.write(f'{clade_}\n')
                res.write(f'\n\n\n')
    result_clades = [clades[1] for clades in suspicious]
    nonredundant_clades = nonredundant(result_clades)
    problematic_orgs, contaminations = problematic(nonredundant_clades)
    plot_problematic(problematic_orgs)
    for tree in trees:
        print(tree)
        tree_to_pdf(tree)


# number of genes
# gene lengths