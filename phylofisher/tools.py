#!/usr/bin/env python
import hashlib
import os
import shutil
from collections import Counter
from datetime import datetime

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd


def parse_aligns(args, input_dir):
    """
    Parses all fasta files in the user-given input directory
    Input: NONE
    Output: dictionary with genes a keys and SeqIO iterator as values
    """
    files = [os.path.join(input_dir, x) for x in os.listdir(input_dir) if x.endswith('.fas')]
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

    my_df = pd.DataFrame(to_df, columns=['UniqueID', 'Orthologs'])
    bin_mat = pd.crosstab(my_df.UniqueID, my_df.Orthologs)

    return bin_mat


def parse_metadata(metadata, input_meta=False):
    """
    Parses metadata.csv to get all org names in each group and subtax
    Input: NONE
    Output: (1) dictionary with groups/subtax as keys and sets of orgs in those groups/subtax as values
            (2) list of all orgs in metadata
    """
    taxa_dict = {}
    with open(metadata, 'r') as infile:
        for line in infile:
            line = line.strip()
            if "Unique ID" not in line:
                if input_meta:
                    _, _, org, group, subtax, _, long_name, _, _ = line.split('\t')
                else:
                    org, long_name, group, subtax, _, _ = line.split('\t')
                taxa_dict[org] = [group, subtax, long_name]

    return taxa_dict


def completeness(args, input_dir, genes=True):
    """
    Computes completeness of genes based on either all taxa_comp_df or a set of taxa_comp_df provided by user
    """

    matrix = parse_aligns(args, input_dir)
    taxa_count, gene_count = matrix.shape
    gene_comp = matrix.sum().divide(other=taxa_count)
    taxa_comp = matrix.sum(axis=1).divide(other=gene_count)

    if genes:
        return matrix
    else:
        return taxa_comp, len(gene_comp)


def make_plot(s, plot_name, y_count, genes=True):
    """
    Creates a plot of gene completeness. Plot is a bar chart with genes sorted from highest to lowest completeness on
        the x-axis and percent complete on the y-axis. Completeness is calculated as percent of taxa_comp_df the gene is present
        in.
    Input: Pandas DataFrame, plot name
    Output: PDF of plot
    """

    if genes:
        calculated_s = 'Ortholog'
        calculated_p = 'Orthologs'
        calc_by = 'Taxa'
    else:
        calculated_s = 'Taxon'
        calculated_p = 'Taxa'
        calc_by = 'Orthologs'

    s = s.sort_values(ascending=False)
    # Dictionary with completeness as key and color as threshold. Color will be used to color bars that are above
    #   threshold.
    colors_threshold = {0.0: 'coral', 0.1: 'black', 0.2: 'lime', 0.3: 'maroon', 0.4: 'khaki',
                        0.5: 'purple', 0.6: 'orange', 0.7: 'olive', 0.8: 'grey', 0.9: 'navy'}
    colors_threshold_inv = {v: k for k, v in colors_threshold.items()}

    # Assigns color to genes depending on the completeness. Based off dict above
    s_colors = [''] * len(s.index)
    for i, row in enumerate(s):
        for j, threshold in enumerate(colors_threshold.keys()):
            if row >= threshold:  # Assigns color based on completeness value
                s_colors[i] = colors_threshold[threshold]
    df = s.to_frame()
    df['colors'] = s_colors  # Adds color list as a column in this DataFrame
    df_colors = [x for x in df['colors']]  # Converts list to tuple because df.plot requires tuple not list
    threshold_counts = Counter(df['colors'])

    # initializing plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    s.plot(width=0.85,
           color=df_colors,
           figsize=(35, 10),
           legend=False,
           fontsize=8,
           ax=ax,
           kind='bar')

    ax.set_title(f'{calculated_s} Completeness by {calc_by}', fontsize=30)

    # Legend
    legend_data = []
    for i, x in enumerate(sorted(colors_threshold.keys(), reverse=True)):
        if i == len(threshold_counts):
            break
        if i == 0:
            label_str = f'Completeness ≥ {round(x * 100)}% ({threshold_counts[colors_threshold[x]]} {calculated_p})'
        else:
            label_str = (f'{round(x * 100)}% ≤ Completeness < {round((x + 0.1) * 100)}%'
                         f' ({threshold_counts[colors_threshold[x]]} {calculated_p})')
        legend_data.append(mpatches.Patch(color=colors_threshold[x],
                                          label=label_str))

    ax.legend(handles=legend_data, facecolor='white', framealpha=1)

    # Y-axis settings and formatting
    ax.set_ylabel(f'Percent Complete ({y_count} {calc_by})', labelpad=15, fontsize=20)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_ylim(0, 1)
    ax.yaxis.set_ticks(np.arange(0, 1.1, 0.1))
    ax.tick_params(axis='y', which='both', labelleft=True, labelright=True)

    # X-axis settingsand formattings
    ax.set_xlabel(f'{calculated_p} ({len(df.index)})', labelpad=15, fontsize=20)

    # Horizontal Lines
    for x in threshold_counts.keys():  # Iterate through colors with counts being the values
        plt.axhline(y=colors_threshold_inv[x], xmin=0, xmax=1, color=x, zorder=0)
    plt.axhline(y=1, xmin=0, xmax=1, color='black', zorder=0)

    # Top X-Axis
    gene_num = len(df.index) - 1
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(range(-1, gene_num, 5))
    ax2.set_xticklabels(range(0, gene_num, 5))
    ax2.xaxis.set_minor_locator(mtick.MultipleLocator(1))

    # Graphic rendering
    fig.set_tight_layout(True)
    fig.savefig(f'{plot_name}.pdf')
    fig.show()


def get_md5(filename):
    # Open,close, read file and calculate MD5 on its contents
    with open(filename, 'rb') as file_to_check:
        # read contents of the file
        data = file_to_check.read()
        # pipe contents of the file through
        md5_returned = hashlib.md5(data).hexdigest()

    return md5_returned


def first_backup(dfo, date_time):
    os.mkdir(f'{dfo}/backups')
    os.mkdir(f'{dfo}/backups/{date_time}')
    shutil.copytree(f'{dfo}/orthologs', f'{dfo}/backups/{date_time}/orthologs')
    shutil.copytree(f'{dfo}/paralogs', f'{dfo}/backups/{date_time}/paralogs')
    shutil.copytree(f'{dfo}/proteomes', f'{dfo}/backups/{date_time}/proteomes')
    shutil.copy(f'{dfo}/metadata.tsv', f'{dfo}/backups/{date_time}/metadata.tsv')
    shutil.copy(f'{dfo}/tree_colors.tsv', f'{dfo}/backups/{date_time}/tree_colors.tsv')


def backup(dfo):
    now = datetime.now().strftime('%d-%b-%Y_%H-%M-%S')

    if os.path.isdir(f'{dfo}/backups'):
        backups = os.listdir(f'{dfo}/backups')
        backups.sort(key=lambda date: datetime.strptime(date, '%d-%b-%Y_%H-%M-%S'), reverse=True)
        latest_backup = backups[0]

        back_files = {}
        db_files = {}
        os.mkdir(f'{dfo}/backups/{now}')
        os.mkdir(f'{dfo}/backups/{now}/orthologs')
        os.mkdir(f'{dfo}/backups/{now}/paralogs')
        os.mkdir(f'{dfo}/backups/{now}/proteomes')

        for root, dirs, files in os.walk(f'{dfo}/backups'):
            for file in files:
                back_files[file] = (os.path.join(root, file))

        for folder in ['orthologs', 'paralogs', 'proteomes']:
            for root, dirs, files in os.walk(f'{dfo}/{folder}'):
                for file in files:
                    db_files[file] = (os.path.join(root, file))

        db_files['metadata.tsv'] = f'{dfo}/metadata.tsv'
        db_files['tree_colors.tsv'] = f'{dfo}/tree_colors.tsv'

        for k, v in db_files.items():
            if k == 'metadata.tsv' or k == 'tree_colors.tsv':
                dest = k
            else:
                dest = '/'.join(v.split('/')[-2:])

            if k in back_files.keys() and get_md5(db_files[k]) == get_md5(back_files[k]):
                os.symlink(f'{dfo}/backups/{latest_backup}/{dest}', f'{dfo}/backups/{now}/{dest}')
            else:
                shutil.copy(v, f'{dfo}/backups/{now}/{dest}')

    else:
        first_backup(dfo, now)
