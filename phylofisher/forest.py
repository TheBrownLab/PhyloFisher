#!/usr/bin/env python
import configparser
import glob
import os
import shutil
import tarfile
import textwrap
from collections import defaultdict, Counter
from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import SeqIO
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
from matplotlib.backends.backend_pdf import PdfPages

from phylofisher import help_formatter

plt.style.use('ggplot')


def configure_colors():
    color_dict = dict()
    with open(color_conf, 'r') as infile:
        infile.readline()
        for line in infile:
            line = line.strip()
            tax, color = line.split('\t')
            color_dict[tax] = color

    return color_dict


def parse_metadata(metadata, input_metadata=None):
    """
    Parse metadata from dataset and input_metadata (if provided)
    input:  metadata csv file, input metadata csv file (optional)
    return: dictionary with combined metadata, dictionary with taxonomy: color
    """
    color_dict = configure_colors()
    metadata_comb = {}
    for line_ in open(metadata):
        if 'Full Name' not in line_:
            sline = line_.split('\t')
            tax = sline[0].strip()
            group = sline[2].strip()
            sub_tax = sline[3]
            full = sline[1].strip()
            if group not in color_dict or color_dict[group].lower() in ['x', 'xx']:
                color_dict[group] = 'white'
            metadata_comb[tax] = {'Higher Taxonomy': group, 'col': color_dict[group], 'full': full,
                                  'Lower Taxonomy': sub_tax}
    if input_metadata:
        for line in open(input_metadata):
            if "FILE_NAME" not in line:
                metadata_input = line.split('\t')
                tax = metadata_input[2].strip().split('_')[0]
                group = metadata_input[3].strip()
                full = metadata_input[6].strip()
                sub_tax = metadata_input[4]
                metadata_comb[tax] = {'Higher Taxonomy': group, 'col': "white", 'full': full, 'Lower Taxonomy': sub_tax}
    return metadata_comb, color_dict


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
            groups.add(metadata[org]['Higher Taxonomy'])
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
            if node.name.count('_') == 3:
                org, __, _, rank = node.name.split('_')
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


def parse_contaminants(file):
    """
    Parse contamination file.
    example:
    ==============================
    Pirisoci	Alveolata	Higher Taxonomy
    DiplDSTH	GoniavonGEN	Unique ID
    Fucucera	Ochrophyta	Lower Taxonomy
    ==============================
    input: tab separated file with contaminants
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
    # check that rank is valid (group, Lower Taxonomy, org)
    assert rank in ['Higher Taxonomy', 'Lower Taxonomy',
                    'Unique ID'], f'{rank} has to be Higher Taxonomy,Lower Taxonomy or org'
    for org in parent.get_leaf_names():
        if (org.count('_') != 4):
            if '..' in org:
                # in case of paralogs
                org = org.split('..')[0]
            else:
                org = org.split('_')[0]
            if rank == 'Higher Taxonomy':
                keywords.add(metadata[org]['Higher Taxonomy'])
            elif rank == 'Lower Taxonomy':
                keywords.add(metadata[org]['Lower Taxonomy'])
            elif rank == 'Unique ID':
                keywords.add(org)
    if keywords:
        if len(keywords) == 1:
            # only in a case that expected neighbourhood is exactly what we
            # are expecting to be. Only one Higher Taxonomy, Lower Taxonomy, tax
            if key in list(keywords)[0]:
                return True
            else:
                return False
        else:
            return False
    return expected_neighborhood(parent.up, cont_key_rank)


def collect_contaminants(tree_file, cont_dict):
    """
    Collect name of all sequences where position on tree corresponds to expected place 
    for a contamination.
    input: tree file, contamination dict from parse_contaminants fucntion
    result: set of proven contaminants, set of proven contamination (same names as in csv result tables)
    """
    t = Tree(tree_file)
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    cont_table_names = set()
    contaminants = set()
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
                        contaminants.add(name)
                        cont_table_names.add(table_name)
        n += 1
    return contaminants, cont_table_names


def get_build_len(name_):
    if os.path.isfile(f'{args.input}/{name_}.trimmed') is True:
        with open(f'{args.input}/{name_}.final', 'r') as build, open(f'{args.input}/{name_}.trimmed', 'r') as trimmed:
            trimmed_len = 0
            len_dict = {}
            for record in SeqIO.parse(trimmed, 'fasta'):
                if len(record.seq) > trimmed_len:
                    trimmed_len = len(record.seq)
                len_dict[record.description] = len(str(record.seq).replace('-', '').replace('X', ''))

            build_len = 0
            for record in SeqIO.parse(build, 'fasta'):
                if len(record.seq) > build_len:
                    build_len = len(record.seq)
        return build_len, len_dict, trimmed_len

    else:
        with open(f'{args.input}/{name_}.final', 'r') as infile:
            build_len = 0
            len_dict = {}
            for record in SeqIO.parse(infile, 'fasta'):
                if len(record.seq) > build_len:
                    build_len = len(record.seq)
                len_dict[record.description] = len(str(record.seq).replace('-', '').replace('X', ''))
        return build_len, len_dict


def tree_to_tsvg(tree_file, contaminants=None, backpropagation=None):
    if contaminants is None:
        contaminants = set()
    tree_base = str(os.path.basename(tree_file))

    # what if they will use somethig different than Raxml? We should make some if statement here maybe.
    name_ = tree_base.split('.')[1]

    if os.path.isfile(f'{args.input}/{name_}.trimmed') is True:
        build_len, len_dict, trimmed_len = get_build_len(name_)
        len_info = f'Final Align Len: {build_len}, Trimmed Align Len: {trimmed_len}'
        len_dict = {k: round(v / trimmed_len, 2) for k, v in len_dict.items()}
    else:
        build_len, len_dict = get_build_len(name_)
        len_info = f'Final Align Len: {build_len}'
        len_dict = {k: round(v / build_len, 2) for k, v in len_dict.items()}

    if not backpropagation:
        table = open(f"{output_folder}/{name_.split('_')[0]}.tsv", 'w')
    else:
        table = open(f"{output_folder}/{name_.split('_')[0]}.tsv", 'r')

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
                format_leaves(backpropagation, contaminants, node, node_style, table, top_ranked, len_dict)
                node.set_style(node_style)

    title_face = TextFace(f'<{name_}  {len_info}, {sus_clades} suspicious clades>', bold=True)
    ts.title.add_face(title_face, column=1)
    t.render(f'{output_folder}/{name_}_tree.svg', tree_style=ts, )
    if not backpropagation:  # what what what?
        table.close()


def format_leaves(backpropagation, contaminants, node, node_style, table, top_ranked, len_dict):
    """This function formats the leaves for the svg file and writes to tsv file"""

    # This parts is about leaves
    empty_face = TextFace("\t" * 20)  # just because of Sefs script
    node.add_face(empty_face, column=2, position="aligned")
    node.add_face(empty_face, column=3, position="aligned")
    org = node.name
    if org.count('_') == 3:  # if org in additions
        org_in_add(backpropagation, contaminants, node, org, table, top_ranked, len_dict)

    elif '..' in org:
        para_from_meta(backpropagation, node, org, table, len_dict)

    else:
        ortho_from_meta(backpropagation, node, node_style, org, table, len_dict)


def org_in_add(backpropagation, contaminants, node, org, table, top_ranked, len_dict):
    """for organisms that are in additions"""
    quality = f'{org.split("_")[-2]}_{org.split("_")[-1]}_{len_dict[org]}'
    org = node.name
    unique_id = org.split('_')[0]
    group = metadata[unique_id]['Higher Taxonomy']
    if group in tax_col:
        tax_name = TextFace(f'[{group} {metadata[unique_id]["Lower Taxonomy"]}]', fgcolor=tax_col[group], bold=True)
    else:
        tax_name = TextFace(f'[{group} {metadata[unique_id]["Lower Taxonomy"]}]', bold=True)
    node.add_face(tax_name, column=1, position="aligned")
    if node.name in contaminants:
        # Seqs to delete
        if not backpropagation:
            table.write(f'{metadata[unique_id]["full"]}_{quality}@{unique_id}\t{group}\td\n')
        delf = TextFace(f'{metadata[unique_id]["full"]}_{quality}@{unique_id}', fgcolor='red')
        node.name = ''
        node.add_face(delf, column=0)
    elif node.name in top_ranked:
        # top ranked guys
        tname = TextFace(f'{metadata[unique_id]["full"]}_{quality}@{unique_id}', bold=True)
        if not backpropagation:
            table.write(f'{metadata[unique_id]["full"]}_{quality}@{unique_id}\t{group}\to\n')
        node.name = ''
        node.add_face(tname, column=0, position='branch-right')
    else:
        if not backpropagation:
            table.write(f'{metadata[unique_id]["full"]}_{quality}@{unique_id}\t{group}\tp\n')
        node.name = f'{metadata[unique_id]["full"]}_{quality}@{unique_id}'


def para_from_meta(backpropagation, node, org, table, len_dict):
    """for paralogs from metadata"""
    unique_id = org.split('..')[0]
    group = f"{metadata[unique_id]['Higher Taxonomy']}"
    paraf = TextFace(f'{metadata[unique_id]["full"]}_{len_dict[org]}@{org}', fgcolor='blue')
    node.name = ''
    node.add_face(paraf, column=0)
    if not backpropagation:
        table.write(f'{metadata[unique_id]["full"]}_{len_dict[org]}@{org}\t{group}\tp\n')
    gface = TextFace(f'[{group} {metadata[unique_id]["Lower Taxonomy"]}]')
    node.add_face(gface, column=1, position="aligned")


def ortho_from_meta(backpropagation, node, node_style, org, table, len_dict):
    """for orthologs from dataset"""
    group = f"{metadata[org]['Higher Taxonomy']}"
    gface = TextFace(f'[{group} {metadata[org]["Lower Taxonomy"]}]')  # TODO do not touch me pleeeease
    color = metadata[org]['col']
    node_style["bgcolor"] = color
    if not backpropagation:
        table.write(f'{metadata[org]["full"]}_{len_dict[org]}@{org}\t{group}\to\n')
    node.name = f'{metadata[org]["full"]}_{len_dict[org]}@{org}'
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
                taxons.add(metadata[org]['Higher Taxonomy'])
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
            taxons.append(metadata[org]["Higher Taxonomy"])
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
            if group_perc[metadata[org]["Higher Taxonomy"]] == 50:
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
            plt.title(f'{org}({metadata[org]["Higher Taxonomy"]})')
            plt.xticks(rotation='vertical')
            plt.yticks(fontsize=4)
            pdf.savefig(bbox_inches='tight')
            plt.close()


def backpropagate_contamination(tree_file, cont_names):
    """
    Changes status in all csv result tables for all proven contaminants
    to 'd' as delete.
    input: tree file, set of proven contaminants (same name format as in csv tables)
    result: None
    """
    tree_base = str(os.path.basename(tree_file))
    output_base = f"{tree_base.split('.')[1].split('_')[0]}"
    tree_base = tree_base.split("_")[0]

    # change original tables before parsing
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

    # changed already parsed table (*_parsed.tsv) if it already exists
    parsed_table_file = f"{output_folder}/{output_base}_parsed.tsv"
    if os.path.isfile(parsed_table_file):
        parsed_table = open(parsed_table_file, 'r').readlines()
        with open(f"{output_folder}/{output_base}.tsv", 'w') as res_:
            for line in parsed_table:
                sline = line.split('\t')
                name = sline[0]
                tax = sline[1]
                status = sline[2].strip()
                if name in cont_names:
                    status = 'd'
                res_.write(f'{name}\t{tax}\t{status}\n')


if __name__ == '__main__':
    description = 'Render .svg and .tsv files of single gene trees for visualization with ParaSorter.'
    parser, optional, required = help_formatter.initialize_argparse(name='forest.py',
                                                                    desc=description,
                                                                    usage='forest.py [OPTIONS] -i <in_dir>')
    # Add Arguments Specific to this script
    # Optional Arguments
    optional.add_argument('-a', '--contaminants', metavar='<contams.tsv>',
                          help=textwrap.dedent("""\
                          Path to file containing previously known contaminants to be removed."""))
    optional.add_argument('-b', '--backpropagate', action='store_true',
                          help=textwrap.dedent("""\
                          Remove contaminants through backpropagation."""))
    optional.add_argument('-t', '--threads', metavar='<N>', type=int,
                          help=textwrap.dedent("""\
                          Number of threads to be used, where N is an integer. 
                          Default: N=1"""))
    optional.add_argument('--local_run', action='store_true',
                          help=textwrap.dedent("""\
                          To be used when sgt_constructor_out.tar.gz has been downloaded from a server for
                          single gene tree visualizations to be rendered locally.
                          """))

    in_help = 'Path to sgt_constructor_out_<M.D.Y>/trees if run remotely or sgt_constructor_out_<M.D.Y>.tar.gz ' \
              'if run locally'
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, in_help=in_help)

    output_folder = args.output

    if args.local_run:
        with tarfile.open(args.input, "r:gz") as tar:
            tar.extractall()
            tar.close()

        args.input = args.input.split('.tar.gz')[0]
        trees_folder = f'{args.input}/trees'
        args.metadata = f'{args.input}/metadata.tsv'
        args.input_metadata = f'{args.input}/input_metadata.tsv'
        color_conf = f'{args.input}/tree_colors.tsv'

        args.input = f'{args.input}/trees'

    else:
        trees_folder = args.input
        config = configparser.ConfigParser()
        config.read('config.ini')
        dfo = str(Path(config['PATHS']['database_folder']).resolve())
        args.metadata = str(os.path.join(dfo, 'metadata.tsv'))
        color_conf = str(Path(config['PATHS']['color_conf']).resolve())
        args.input_metadata = str(os.path.abspath(config['PATHS']['input_file']))

    if not args.backpropagate:
        os.mkdir(output_folder)

    trees = glob.glob(f"{trees_folder}/*.tre")

    number_of_genes = len(trees)
    metadata, tax_col = parse_metadata(args.metadata, args.input_metadata)
    threads = args.threads

    if not args.backpropagate:
        suspicious = parallel_susp_clades(trees)
        suspicious = sorted(suspicious)

        make_txt = False
        for x, y in suspicious:
            if y:
                make_txt = True

        if make_txt:
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

        # Commented out due to issues with the function
        # plot_major_taxa(orgs_taxa)

    if args.contaminants:
        cont_dict = parse_contaminants(args.contaminants)
        for tree in trees:
            contaminants, contaminated_table = collect_contaminants(tree, cont_dict)
            if args.backpropagate:
                backpropagate_contamination(tree, contaminated_table)
                tree_to_tsvg(tree, contaminants, backpropagation=True)
            else:
                tree_to_tsvg(tree, contaminants)
    else:
        for tree in trees:
            tree_to_tsvg(tree)

    if args.local_run:
        shutil.rmtree('/'.join(args.input.split('/')[:-1]))
