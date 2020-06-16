#!/usr/bin/env python
import argparse
import os
import shutil
import sys
import textwrap
from collections import defaultdict
from pathlib import Path
import configparser

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patches as mpatches
import numpy as np
from phylofisher import help_formatter


def parse_metadata():
    """
    Parses metadata.csv to get all org names in each group and subtax
    Input: NONE
    Output: (1) dictionary with groups/subtax as keys and sets of orgs in those groups/subtax as values
            (2) list of all orgs in metadata
    """
    groups = defaultdict(set)
    all_orgs = []
    for line in open(metadata):
        if line.strip():
            if "Full Name" not in line:
                org, _, group, subtax, _, _ = line.split('\t')
                groups[group].add(org)
                groups[subtax].add(org)
                all_orgs.append(org)

    return groups, all_orgs


def parse_taxa_list():
    """
    Parses user-input taxa list
    Input: NONE
    Ouput: List of user-input taxa
    """
    in_taxa = []
    possible_orgs = set()
    groups, _ = parse_metadata()
    with open(args.taxa, 'r') as infile:
        for line in infile:
            in_taxa.append(line.strip())

    for taxon in in_taxa:
        if taxon in groups.keys():
            orgs = [org for org in groups[taxon]]
            for x in orgs:
                possible_orgs.add(x)
        else:
            possible_orgs.add(taxon)

    # Needs to return all possible orgs
    return possible_orgs


def parse_aligns():
    """
    Parses all fasta files (w/ given suffix if provided) in the user-given input directory
    Input: NONE
    Output: dictionary with genes a keys and SeqIO iterator as values
    """
    files = [os.path.join(args.input, x) for x in os.listdir(args.input) if x.endswith(args.suffix)]
    all_taxa_records = {}
    for file in files:
        gene = os.path.basename(file).split('.')[0]
        with open(file, 'r') as infile:
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    taxon = line[1:]
                    if len(taxon.split('_')) == 4:
                        taxon = taxon.split('_')[0]

                    if taxon not in all_taxa_records.keys():
                        all_taxa_records[taxon] = [gene]
                    else:
                        all_taxa_records[taxon].append(gene)

    to_df = []
    for taxon in all_taxa_records:
        for gene in all_taxa_records[taxon]:
            to_df.append([taxon, gene])

    my_df = pd.DataFrame(to_df, columns=['Taxa', 'Gene'])
    bin_mat = pd.crosstab(my_df.Taxa, my_df.Gene)

    return bin_mat


def completeness():
    """
    Comptutes completeness of genes based on either all taxa or a set of taxa provided by user
    """
    matrix = parse_aligns()
    if args.taxa:
        possible_orgs = parse_taxa_list()
    else:
        _, possible_orgs = parse_metadata()

    taxa_count, gene_count = matrix.shape
    gene_comp = matrix.sum().divide(other=taxa_count)
    taxa_comp = matrix.sum(axis=1).divide(other=gene_count)

    return gene_comp, taxa_comp


def make_gene_plot(df, plot_name, taxa_count):
    """
    Creates a plot of gene completeness. Plot is a bar chart with genes sorted from highest to lowest completeness on
        the x-axis and percent complete on the y-axis. Completeness is calculated as percent of taxa the gene is present
        in.
    Input: Pandas DataFrame, plot name
    Output: PDF of plot
    """
    df = df.sort_values(by="% complete", ascending=False)
    # Dictionary with completeness as key and color as threshold. Color will be used to color bars that are above
    #   threshold.
    colors_threshold = {0.0: 'coral', 0.1: 'azure', 0.2: 'lime', 0.3: 'maroon', 0.4: 'khaki',
                        0.5: 'purple', 0.6: 'orange', 0.7: 'olive', 0.8: 'grey', 0.9: 'navy'}
    colors_threshold_inv = {v: k for k, v in colors_threshold.items()}

    # Assigns color to genes depending on the completeness. Based off dict above
    df_colors = [''] * len(df.index)
    for i, row in enumerate(df["% complete"]):
        for j, threshhold in enumerate(colors_threshold.keys()):
            if row >= threshhold:  # Assisns color based on completeness value
                df_colors[i] = colors_threshold[threshhold]
    df['colors'] = df_colors  # Adds color list as a column in this DataFrame
    df_colors = tuple([x for x in df['colors']])  # Converts list to tuple because df.plot requires tuple not list
    threshold_counts = df['colors'].value_counts().to_dict()

    # initailizng plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    df.plot(kind='bar',
            width=0.85,
            color=[df_colors],
            figsize=(35, 10),
            legend=False,
            fontsize=8,
            ax=ax)
    ax.set_title('Gene Completeness by Taxa', fontsize=30)

    # Legend
    legend_data = []
    for i, x in enumerate(sorted(colors_threshold.keys(), reverse=True)):
        if i == len(threshold_counts):
            break
        if i == 0:
            label_str = f'Completeness ≥ {round(x * 100)}% ({threshold_counts[colors_threshold[x]]} Genes)'
        else:
            label_str = (f'{round(x * 100)}% ≤ Completeness < {round((x + 0.1) * 100)}%'
                         f' ({threshold_counts[colors_threshold[x]]} Genes)')
        legend_data.append(mpatches.Patch(color=colors_threshold[x],
                                          label=label_str))

    ax.legend(handles=legend_data, facecolor='white', framealpha=1)

    # Y-axis settings and formatting
    ax.set_ylabel(f'Percent Complete ({taxa_count} taxa)', labelpad=15, fontsize=20)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(start, end + 0.1, 0.1))
    ax.tick_params(axis='y', which='both', labelleft=True, labelright=True)

    # X-axis settingsand formattings
    ax.set_xlabel(f'Genes ({len(df.index)})', labelpad=15, fontsize=20)

    # Horizontal Lines
    for x in threshold_counts.keys():  # Iterate through colors with counts being the values
        plt.axhline(y=colors_threshold_inv[x], xmin=0, xmax=1, color=x, zorder=0)
    plt.axhline(y=1, xmin=0, xmax=1, color='black', zorder=0)

    # Top X-Axis
    gene_num = len(df.index)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(range(-1, gene_num, 5))
    ax2.set_xticklabels(range(0, gene_num, 5))
    ax2.xaxis.set_minor_locator(mtick.MultipleLocator(1))

    # Graphic rendering
    fig.set_tight_layout(True)
    fig.savefig(f'{plot_name}.pdf')
    fig.show()


def genes_to_keep(df):
    genes = []
    df = df.sort_values(["% complete"], ascending=False)
    if args.gene_number:
        count = 0
        for index, _ in df.iterrows():
            if count < args.gene_number:
                genes.append(index)
                count += 1
    elif args.percent_complete:
        for index, row in df.iterrows():
            if row['% complete'] >= args.percent_complete:
                genes.append(index)

    return genes


def subsetter(df):
    if os.path.exists(subset_dir):
        shutil.rmtree(subset_dir)
    os.makedirs(subset_dir)

    genes = genes_to_keep(df)
    files = [os.path.join(args.input, x) for x in os.listdir(args.input) if x.endswith(args.suffix)]
    for gene in genes:
        for file in files:
            if gene == os.path.basename(file).split('.')[0]:
                src = file
                dest = f'{subset_dir}/{os.path.basename(file)}'
                shutil.copy(src, dest)


if __name__ == '__main__':
    description = 'Subsets gene based on gene completeness.'
    parser, optional, required = help_formatter.initialize_argparse(name='missing_data.py',
                                                                    desc=description,
                                                                    usage='missing_data.py '
                                                                          '[OPTIONS] -i <input>')

    # Optional Arguments
    optional.add_argument('-t', '--taxa', metavar='<taxa.tsv>',
                          help=textwrap.dedent("""\
                          TSV file of taxa to be considered for completeness.
                          Default is all taxa
                          """))
    optional.add_argument('--subset', action='store_true',
                          help=textwrap.dedent("""\
                          Subset input genes
                          """))
    optional.add_argument('-n', '--gene_number', type=int, metavar='<N>',
                          help=textwrap.dedent("""\
                          Number of genes in subset.
                          This will be ignored if not used with --subset.
                          Cannot be used with percent_complete.
                          """))
    optional.add_argument('-c', '--percent_complete', type=float, metavar='<N>',
                          help=textwrap.dedent("""\
                          Threshold for percent missing when subsetting.
                          This will be ignored if not used with --subset.
                          Cannot be used with gene_number.
                          """))

    args = help_formatter.get_args(parser, optional, required)

    if args.subset and args.gene_number and args.percent_complete:
        parser.error("--subset requires either --gene_number or --percent_complete NOT both.")
    elif args.subset and (not args.gene_number or not args.percent_complete):
        parser.error("--subset requires either --gene_number or --percent_complete.")

    config = configparser.ConfigParser()
    config.read('config.ini')
    dfo = str(Path(config['PATHS']['dataset_folder']).resolve())
    metadata = f'{dfo}/metadata.tsv'

    if os.path.isdir(args.output) is False:
        os.mkdir(args.output)

    completeness_df, taxa = completeness()

    # plot_name = f'{args.output}/out'
    #
    # make_gene_plot(completeness_df, plot_name, len(taxa))
    #
    # if args.subset:
    #     if args.gene_number:
    #         subset_dir = f'{args.output}/subset_n{args.gene_number}'
    #     else:
    #         subset_dir = f'{args.output}/subset_c{args.percent_complete}'
    #     subsetter(df=completeness_df)
