#!/usr/bin/env python
import glob
import os
import textwrap
from collections import defaultdict, Counter
import configparser
from multiprocessing import Pool

from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from phylofisher import help_formatter

plt.style.use('ggplot')


def parse_metadata(metadata, input_metadata=None):
    """
    Parse metadata from dataset and input_metadata (if provided)
    input:  metadata csv file, input metadata csv file (optional)
    return: dictionary with combined metadata, dictionary with taxonomy: color
    """
    metadata_comb = {}
    tax_col = {}
    for line_ in open(metadata):
        if 'Full Name' not in line_:
            sline = line_.split('\t')
            tax = sline[0].strip()
            group = sline[2].strip()
            col = sline[4].strip()
            sub_tax = sline[3]
            full = sline[1].strip()
            if group not in tax_col:
                if col.lower() in ['x', 'xx']:
                    col = 'white'
                tax_col[group] = col
            metadata_comb[tax] = {'group': group, 'col': tax_col[group], 'full': full, 'subtax': sub_tax}
    if input_metadata:
        for line in open(input_metadata):
            if "FILE_NAME" not in line:
                metadata_input = line.split('\t')
                tax = metadata_input[2].strip().split('_')[0]
                group = metadata_input[3].strip()
                full = metadata_input[6].strip()
                sub_tax = metadata_input[4]
                metadata_comb[tax] = {'group': group, 'col': "white", 'full': full, 'subtax': sub_tax}
    return metadata_comb, tax_col


def suspicious_clades(tree):
    """
    Find suspicious clades (more than 70 bs and more than 2 tax groups)
    input: phylogenetic tree
    output: tuple of tree name and list of suspicious clades 
    """
    t = Tree(tree)
    # midpoint rooted tree
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)

    supported_clades = []
    for node in t.traverse('preorder'):
        if (node.is_root() is False) and (node.is_leaf() is False):
            # report only clades which encompass less than a half of all oranisms
            if node.support >= 70 and (len(node) < (len(t) - len(node))):
                clade = node.get_leaf_names()
                if len(clade) > 1:  # do we need this statement?
                    supported_clades.append(clade)
    suspicious = []
    for clade in supported_clades:
        groups = set()
        for org in clade:
            # get org name
            if '..' in org:
                org = org.split('..')[0]
            else:
                org = org.split('_')[0]
            groups.add(metadata[org]['group'])
        if len(groups) > 1:
            suspicious.append(clade)
    return tree, suspicious


def get_best_candidates(tree_file):
    """
    Check best candidates after trimming.
    example: if q1 doesn't survived -> q2 is ranked as best candidate == ortholog
    input: tree file
    output: set of best candidate sequences
    """
    t = Tree(tree_file)
    top_rank = defaultdict(dict)
    for node in t.traverse('preorder'):
        if node.is_leaf():
            if node.name.count('_') == 4:
                org, __, _, rank, _ = node.name.split('_')
                rank = int(rank[1:-1])
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


def parse_contaminations(file):
    """
    Parse contamination file.
    example:
    ==============================
    Pirisoci	Alveolata	group
    DiplDSTH	GoniavonGEN	org
    Fucucera	Ochrophyta	subtax
    ==============================
    input: tab separated file with contaminations
    result: dictionary with organism as key and tuple with taxonomy
            and rank 
    """
    cont_dict = {}
    for line in open(file):
        org, tax, rank = line.split('\t')
        cont_dict[org] = (tax, rank.strip())
    return cont_dict


def expected_neighborhood(parent, cont_key_rank):
    """
    Recursively check if ornanisms around target has some taxonomical group 
    (not added in this run of fisher) and evaluete if target has expected 
    neighbourhood.
    input: parent node, (key, rank) => example: (Alveolata, group) or (Gonianon, org)
    output: True/False
    """
    keywords = set()
    key = cont_key_rank[0]
    rank = cont_key_rank[1]
    # check that rank is valid (group, subtax, org)
    assert rank in ['group', 'subtax', 'org'], f'{rank} has to be group,subtax or org'
    for org in parent.get_leaf_names():
        if (org.count('_') != 4):
            if '..' in org:
                # in case of paralogs
                org = org.split('..')[0]
            else:
                org = org.split('_')[0]
            if rank == 'group':
                keywords.add(metadata[org]['group'])
            elif rank == 'subtax':
                keywords.add(metadata[org]['subtax'])
            elif rank == 'org':
                keywords.add(org)
    if keywords:
        if len(keywords) == 1:
            # only in a case that expected neighbourhood is exactly what we
            # are expecting to be. Only one group, subtax, tax
            if key in list(keywords)[0]:
                return True
            else:
                return False
        else:
            return False
    return expected_neighborhood(parent.up, cont_key_rank)


def collect_contaminations(tree_file, cont_dict):
    """
    Collect name of all sequences where position on tree corresponds to expected place 
    for a contamination.
    input: tree file, contamination dict from parse_contaminations fucntion
    result: set of proven contaminations, set of proven contamination (same names as in csv result tables)
    """
    t = Tree(tree_file)
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    cont_table_names = set()
    contaminations = set()
    n = 0
    for node in t.traverse('preorder'):
        if node.is_leaf() is True:
            if node.name.count('_') == 4:
                name = node.name
                org = name.split('_')[0]
                quality = f'{node.name.split("_")[-3]}_{node.name.split("_")[-2]}_{node.name.split("_")[-1]}'
                table_name = f'{metadata[org]["full"]}_{quality}@{org}'
                if org in cont_dict:
                    exp_hood = expected_neighborhood(node.up, cont_dict[org])
                    if exp_hood is True:
                        contaminations.add(name)
                        cont_table_names.add(table_name)
        n += 1
    return contaminations, cont_table_names


def tree_to_tsvg(tree_file, contaminations=None, backpropagation=None):
    if contaminations is None:
        contaminations = set()
    tree_base = str(os.path.basename(tree_file))
    if args.prefix:
        tree_base = tree_base.replace(args.prefix, '')
    if args.suffix:
        tree_base = tree_base.replace(args.suffix, '')

    output_base = f"{tree_base.split('.')[1].split('_')[0]}"
    if not backpropagation:
        table = open(f"{output_folder}/{output_base.split('_')[0]}.tsv", 'w')
    else:
        table = open(f"{output_folder}/{output_base.split('_')[0]}.tsv", 'r')

    top_ranked = get_best_candidates(tree_file)
    t = Tree(tree_file)
    ts = TreeStyle()
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    sus_clades = 0

    for node in t.traverse('preorder'):
        node_style = NodeStyle()
        node_style['vt_line_width'] = 3
        node_style['hz_line_width'] = 3
        node_style['vt_line_type'] = 0
        node_style['hz_line_type'] = 0

        if node.is_root() is False:
            if node.is_leaf() is False:
                # All internal nodes
                supp, sus_clades = format_nodes(node, node_style, sus_clades, t)
                node.add_face(supp, column=0, position="branch-bottom")
                node.set_style(node_style)
            else:
                # All leaves
                format_leaves(backpropagation, contaminations, node, node_style, table, top_ranked)
                node.set_style(node_style)

    name_, trim_len = tree_base.split('.')[1].split('_')
    title_face = TextFace(f'<{name_}  trim_aln_len: {trim_len}, {sus_clades} suspicious clades>', bold=True)
    ts.title.add_face(title_face, column=1)
    t.render(output_base + '_tree.svg', tree_style=ts)
    if not backpropagation:  # what what what?
        table.close()


def format_leaves(backpropagation, contaminations, node, node_style, table, top_ranked):
    """This function formats the leaves for the svg file and writes to tsv file"""

    # This parts is about leaves
    empty_face = TextFace("\t" * 20)  # just because of Sefs script
    node.add_face(empty_face, column=2, position="aligned")
    node.add_face(empty_face, column=3, position="aligned")
    org = node.name
    if org.count('_') == 4:  # if org in additions
        org_in_add(backpropagation, contaminations, node, org, table, top_ranked)

    elif '..' in org:
        para_from_meta(backpropagation, node, org, table)

    else:
        ortho_from_meta(backpropagation, node, node_style, org, table)


def org_in_add(backpropagation, contaminations, node, org, table, top_ranked):
    """for organisms that are in additions"""
    quality = f'{org.split("_")[-3]}_{org.split("_")[-2]}_{org.split("_")[-1]}'
    org = node.name.split('_')[0]
    group = metadata[org]['group']
    if group in tax_col:
        tax_name = TextFace(f'[{group} {metadata[org]["subtax"]}]', fgcolor=tax_col[group], bold=True)
    else:
        tax_name = TextFace(f'[{group} {metadata[org]["subtax"]}]', bold=True)
    node.add_face(tax_name, column=1, position="aligned")
    if node.name in contaminations:
        # Seqs to delete
        if not backpropagation:
            table.write(f'{metadata[org]["full"]}_{quality}@{org}\t{group}\td\n')
        delf = TextFace(f'{metadata[org]["full"]}_{quality}@{org}', fgcolor='red')
        node.name = ''
        node.add_face(delf, column=0)
    elif node.name in top_ranked:
        # top ranked guys
        tname = TextFace(f'{metadata[org]["full"]}_{quality}@{org}', bold=True)
        if not backpropagation:
            table.write(f'{metadata[org]["full"]}_{quality}@{org}\t{group}\to\n')
        node.name = ''
        node.add_face(tname, column=0, position='branch-right')
    else:
        if not backpropagation:
            table.write(f'{metadata[org]["full"]}_{quality}@{org}\t{group}\tp\n')
        node.name = f'{metadata[org]["full"]}_{quality}@{org}'


def para_from_meta(backpropagation, node, org, table):
    """for paralogs from metadata"""
    para, length = org.split('_')
    org = para.split('..')[0]
    group = f"{metadata[org]['group']}"
    paraf = TextFace(f'{metadata[org]["full"]}_{length}@{para}', fgcolor='blue')
    node.name = ''
    node.add_face(paraf, column=0)
    if not backpropagation:
        table.write(f'{metadata[org]["full"]}_{length}@{para}\t{group}\tp\n')
    gface = TextFace(f'[{group} {metadata[org]["subtax"]}]')
    node.add_face(gface, column=1, position="aligned")


def ortho_from_meta(backpropagation, node, node_style, org, table):
    """for orthologs from dataset"""
    org, length = org.split('_')
    group = f"{metadata[org]['group']}"
    gface = TextFace(f'[{group} {metadata[org]["subtax"]}]')  # TODO do not touch me pleeeease
    color = metadata[org]['col']
    node_style["bgcolor"] = color
    if not backpropagation:
        table.write(f'{metadata[org]["full"]}_{length}@{org}\t{group}\to\n')
    node.name = f'{metadata[org]["full"]}_{length}@{org}'
    node.add_face(gface, column=1, position="aligned")


def format_nodes(node, node_style, sus_clades, t):
    """This function visually formats the nodes in the svg file based on the nodes' support values and whether or not
    the clade is suspicous"""
    supp = TextFace(f'{int(node.support)}', fsize=8)
    if node.support >= 70:
        supp.bold = True
        taxons = set()
        orgs = node.get_leaf_names()
        if len(orgs) > 1:
            for org in orgs:
                if '..' in org:  # for paralogs
                    org = org.split('..')[0]
                else:
                    org = org.split('_')[0]  # for potential orthologs
                taxons.add(metadata[org]['group'])
        if len(taxons) > 1 and (len(node) < (len(t) / 2)):
            node_style['shape'] = 'sphere'
            node_style['size'] = 12
            node_style['fgcolor'] = 'red'
            node_style['bgcolor'] = 'Silver'
            sus_clades += 1
    else:
        supp.fsize = 7
    return supp, sus_clades


def parallel_susp_clades(trees):
    """Parallelizes the function suspicious_clades()"""
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
                    for other in sorted_sets[i + 1:]:
                        if set_.issubset(other):
                            redundant = True
                if redundant is False:
                    nonredundant.append(set_)
    return nonredundant


def collect_major_taxa(nonredundant_clades_):
    """
    Collect information about major taxonomic group for all nonreduntant clades
    and connect this information with all organisms in that clade. This function
    prepares data for plot_major_taxa function.
    input: set of nonredundant clades
    return: dictionary with orgs as keys and list of taxonomic groups (major tax groups
    from clades with a given organism) as values
    """
    orgs_ = defaultdict(list)
    for clade in nonredundant_clades_:
        # collect all taxonomic groups for a given clade
        taxons = []
        for org in clade:
            if '..' in org:
                # for paralogs
                org = org.split('..')[0]
            else:
                org = org.split('_')[0]
            taxons.append(metadata[org]["group"])
        group_perc = {}
        for tax in set(taxons):
            group_perc[tax] = (taxons.count(tax) / len(taxons)) * 100
        max_perc = max(group_perc, key=group_perc.get)
        for seq in clade:
            # parse names again
            if '..' in seq:
                # paralogs
                org = seq.split('..')[0]
            else:
                org = seq.split('_')[0]
            if group_perc[metadata[org]["group"]] == 50:
                orgs_[org].append('unspecified')
            else:
                orgs_[org].append(max_perc)
    return orgs_


def plot_major_taxa(orgs_taxa):
    """
    Plot major taxonomic groups from clades (with a given org) for all organisms.
    input: result from collect_major_taxa function (dictionary with list as keys)
    return: None
    """
    max_ = 0
    for org, groups in sorted(orgs_taxa.items()):
        # this 'for loop' is just looking for maximum value of clades
        # with same taxonomic group for one organism. This number (max_)
        # is then used as maximum in plt.ylim to put all orgs on the same
        # scale
        counted_groups = dict(Counter(groups))
        group_total = 0
        for group in counted_groups.values():
            group_total += group
        if group_total > max_:
            max_ = group_total
    with PdfPages('orgs_taxa.pdf') as pdf:
        for org, groups in sorted(orgs_taxa.items()):
            counted_groups = dict(Counter(groups))
            colors = []
            for group in counted_groups.keys():
                if group == 'unspecified':
                    colors.append('black')
                else:
                    colors.append(tax_col.get(group, "white"))
            plt.bar(counted_groups.keys(), counted_groups.values(), color=colors)
            plt.tight_layout()
            plt.ylim(0, max_)
            plt.title(f'{org}({metadata[org]["group"]})')
            plt.xticks(rotation='vertical')
            plt.yticks(fontsize=4)
            pdf.savefig(bbox_inches='tight')
            plt.close()


def backpropagate_contamination(tree_file, cont_names):
    """
    Changes status in all csv result tables for all proven contaminations
    to 'd' as delete.
    input: tree file, set of proven contaminations (same name format as in csv tables)
    result: None
    """
    tree_base = str(os.path.basename(tree_file))
    output_base = f"{tree_base.split('.')[1].split('_')[0]}"
    if args.prefix:
        tree_base = tree_base.replace(args.prefix, '')
    if args.suffix:
        tree_base = tree_base.replace(args.suffix, '')
    tree_base = tree_base.split("_")[0]
    original_table = open(f"{output_folder}/{output_base}.tsv", 'r').readlines()
    with open(f"{output_folder}/{output_base}.tsv", 'w') as res_:
        for line in original_table:
            sline = line.split('\t')
            name = sline[0]
            tax = sline[1]
            status = sline[2].strip()
            if name in cont_names:
                status = 'd'
            res_.write(f'{name}\t{tax}\t{status}\n')


if __name__ == '__main__':
    description = 'Inspects single gene trees for contamination.'
    parser, optional, required = help_formatter.initialize_argparse(name='forest.py',
                                                                    desc=description,
                                                                    usage='forest.py [OPTIONS] -i <in_dir>',
                                                                    dataset=True,
                                                                    input_meta=True)
    # Add Arguments Specific to this script
    # Optional Arguments
    optional.add_argument('-a', '--contaminations', metavar='<contams>',
                          help=textwrap.dedent("""\
                          Path to file containing previously known contaminations to be removed."""))
    optional.add_argument('-b', '--backpropagate', action='store_true',
                          help=textwrap.dedent("""\
                          Remove contaminations through backpropagation."""))
    optional.add_argument('-t', '--threads', metavar='<N>',
                          help=textwrap.dedent("""\
                          Number of threads to be used, where N is an integer. 
                          Default: N=1"""))

    args = help_formatter.get_args(parser, optional, required)

    trees_folder = args.input
    output_folder = args.output

    if args.dataset_folder:
        dfo = str(Path(args.dataset_folder))
    else:
        config = configparser.ConfigParser()
        config.read('config.ini')
        dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
    args.metadata = str(os.path.join(dfo, 'metadata.tsv'))

    if not args.input_metadata:
        config = configparser.ConfigParser()
        config.read('config.ini')
        args.input_metadata = str(os.path.abspath(config['PATHS']['input_file']))

    if not args.backpropagate:
        os.mkdir(output_folder)

    trees = glob.glob(f"{trees_folder}/{args.prefix}*{args.suffix}")

    number_of_genes = len(trees)
    metadata, tax_col = parse_metadata(args.metadata, args.input_metadata)
    threads = args.threads

    if not args.backpropagate:
        suspicious = parallel_susp_clades(trees)
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
        orgs_taxa = collect_major_taxa(nonredundant_clades)
        plot_major_taxa(orgs_taxa)

    if args.contaminations:
        cont_dict = parse_contaminations(args.contaminations)
        for tree in trees:
            contaminations, contaminated_table = collect_contaminations(tree, cont_dict)
            if args.backpropagate:
                backpropagate_contamination(tree, contaminated_table)
                tree_to_tsvg(tree, contaminations, backpropagation=True)
            else:
                tree_to_tsvg(tree, contaminations)
    else:
        for tree in trees:
            tree_to_tsvg(tree)
