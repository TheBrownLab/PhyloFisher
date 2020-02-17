#!/usr/bin/env python
import os
import argparse
import glob
import shutil
import textwrap
from collections import defaultdict
from Bio import SeqIO
import statistics
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
from shutil import copyfile
import configparser
from phylofisher import fisher


def parse_metadata():
    """
    Parses metadata.csv to get all org names in each group and subtax
    Input: NONE
    Output: (1)dictionary with groups/subtax as keys and sets of orgs in those groups/subtax as values
            (2) list of all orgs in metadata
    """
    groups = defaultdict(set)
    all_orgs = []
    for line in open(args.metadata):
        if line.strip():
            if "Full Name" not in line:
                org, _, group, subtax, _, _, _ = line.split('\t')
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
    all_gene_records = {}
    for file in files:
        gene = os.path.basename(file).split('.')[0]
        records = SeqIO.parse(open(file, 'r'), 'fasta')
        all_gene_records[gene] = records

    return all_gene_records


def completeness():
    """
    Comptutes completeness of genes based on either all taxa or a seet of taxa provided by user
    """
    if args.taxa:
        possible_orgs = parse_taxa_list()
    else:
        _, possible_orgs = parse_metadata()

    gene_count = defaultdict(int)
    gene_comp = defaultdict(float)
    gene_records = parse_aligns()
    my_orgs = set()

    for gene in gene_records:
        for record in gene_records[gene]:
            if record.id in possible_orgs:
                gene_count[gene] += 1
                my_orgs.add(record.id)

    for gene in gene_count:
        gene_comp[gene] = gene_count[gene] / len(my_orgs)

    df = pd.DataFrame.from_dict(gene_comp, orient='index', columns=["% complete"])
    return df, my_orgs


def make_plot(df, plot_name, taxa_count):
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
            if row >= threshhold:                                   # Assisns color based on completeness value
                df_colors[i] = colors_threshold[threshhold]
    df['colors'] = df_colors                                        # Adds color list as a column in this DataFrame
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
        else :
            label_str = f'{round(x * 100)}% ≤ Completeness < {round((x + 0.1) * 100)}% ({threshold_counts[colors_threshold[x]]} Genes)'
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
    if os.path.exists(args.output):
        shutil.rmtree(args.output)
    os.makedirs(args.output)

    genes = genes_to_keep(df)
    files = [os.path.join(args.input, x) for x in os.listdir(args.input) if x.endswith(args.suffix)]
    for gene in genes:
        for file in files:
            if gene == os.path.basename(file).split('.')[0]:
                src = file
                dest = f'{args.output}/{os.path.basename(file)}'
                shutil.copy(src, dest)


if __name__ == '__main__':
    formatter = lambda prog: fisher.myHelpFormatter(prog, max_help_position=100)

    parser = argparse.ArgumentParser(prog='missing_data.py',
                                     # TODO: Description
                                     description='some description',
                                     usage='missing_data.py [OPTIONS] -i <in_dir> -t <taxa> '
                                           '-d <dataset> -n <gene_number>',
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent("""\
                                         additional information:
                                            
                                            """))
    optional = parser._action_groups.pop()
    mut_excl = parser.add_mutually_exclusive_group()
    required = parser.add_argument_group('required arguments')

    # Required Arguments
    required.add_argument('-i', '--input', required=True, type=str, metavar='<in_dir>',
                          help=textwrap.dedent("""\
                                  Path to input directory containing gene alignments.
                                  """))

    required.add_argument('-m', '--metadata', required=True, metavar='<dataset>',
                          help=textwrap.dedent("""\
                                  Path to metadata.tsv
                                  """))
    # Mutually Exclusive Arguments
    mut_excl.add_argument('-n', '--gene_number', type=int, metavar='<N>',
                          help=textwrap.dedent("""\
                                  Number of genes for analysis
                                  """))
    mut_excl.add_argument('-c', '--percent_complete', type=float, metavar='<N>',
                          help=textwrap.dedent("""\
                                  Threshold for percent missing
                                  """))
    # Optional Arguments
    optional.add_argument('-o', '--output', type=str, default='missing_data_out', metavar='<out_dir>',
                          help=textwrap.dedent("""\
                                  Path to output directory
                                  """))
    optional.add_argument('-t', '--taxa', metavar='<taxa.tsv>',
                          help=textwrap.dedent("""\
                                  TSV file of taxa to be considered for completeness.
                                  Default is all taxa
                                  """))
    optional.add_argument('-s', '--suffix', type=str, default='', metavar='<suff>',
                          help=textwrap.dedent("""\
                                  Suffix of input files
                                  Default: NONE
                                  Example: path/to/input/*.suffix
                                  """))
    optional.add_argument('-p', '--plot', metavar='taxa_completeness.pdf',
                          help=textwrap.dedent("""\
                                  Plot missing data statistics.
                                  Default: out.pdf
                                  """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                                  Show this help message and exit.
                                  """))

    parser._action_groups.append(optional)
    args = parser.parse_args()

    completeness_df, taxa = completeness()
    if args.plot:
        make_plot(completeness_df, args.plot, len(taxa))

    subsetter(df=completeness_df)
