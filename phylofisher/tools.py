#!/usr/bin/env python
import hashlib
import os
import shutil
from collections import Counter
from datetime import datetime
from peewee import *
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
from phylofisher.db_map import database, BaseModel, Genes, Taxonomies, Metadata, Sequences


def parse_aligns(ortholog=True):
    '''
    Parses the Sequences table in the database to create a binary matrix of taxa and orthologs.

    :param ortholog: whether to include only orthologs (True) or paralogs (False), defaults to True
    :type ortholog: bool, optional
    :return: a binary matrix of completeness
    :rtype: pd.DataFrame
    '''
    
    all_taxa_records = {}

    if ortholog:
        db_query = Sequences.select(Sequences.header, Sequences.gene, Sequences.metadata).where(Sequences.is_paralog == False)
    else:
        db_query = Sequences.select(Sequences.header, Sequences.gene, Sequences.metadata).where(Sequences.is_paralog == True)

    for q in db_query:
        taxon = Metadata.get(Metadata.id == q.metadata).short_name
        gene = Genes.get(Genes.id == q.gene).name

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


def parse_metadata(meta_path, input_meta=False):
    '''
    Parses the input_metadata or database to create a dictionary of taxa and their taxonomic information.

    :param meta_path: path to the metadata file or database
    :type meta_path: str
    :param input_meta: whether the input is a metadata file (True) or database (False), defaults to False
    :type input_meta: bool, optional
    :return: a dictionary of taxa and their taxonomic information
    :rtype: dict
    '''
    taxa_dict = {}
    if input_meta:
        with open(meta_path, 'r') as infile:
            for line in infile:
                line = line.strip()
                if "Unique ID" not in line:
                    _, _, org, group, subtax, _, long_name, _, _ = line.split('\t')
                    taxa_dict[org] = [group, subtax, long_name]
    else:
        database.init(meta_path)
        database.connect()
        db_query = Metadata.select(Metadata.short_name, Metadata.higher_taxonomy, Metadata.lower_taxonomy, Metadata.long_name)
        for q in db_query:
            org = q.short_name
            long_name = q.long_name
            group = Taxonomies.get(Taxonomies.id == q.higher_taxonomy).taxonomy
            subtax = Taxonomies.get(Taxonomies.id == q.lower_taxonomy).taxonomy
            taxa_dict[org] = [group, subtax, long_name]
        database.close()

    return taxa_dict


def completeness(orthologs=True, genes=True):
    '''
    Calculates the completeness of genes or taxa based on the presence of orthologs in the database.

    :param orthologs: whether to calculate completeness based on orthologs (True) or paralogs (False), defaults to True
    :type orthologs: bool, optional
    :param genes: whether to return the completeness of genes (True) or taxa (False), defaults to True
    :type genes: bool, optional
    :return: a tuple containing either a DataFrame of gene completeness or a Series of taxa completeness and the number of genes
    :rtype: tuple
    '''

    matrix = parse_aligns(orthologs)
    taxa_count, gene_count = matrix.shape
    gene_comp = matrix.sum().divide(other=taxa_count)
    taxa_comp = matrix.sum(axis=1).divide(other=gene_count)

    if genes:
        return matrix
    else:
        return taxa_comp, len(gene_comp)


def make_plot(s, plot_name, y_count, genes=True):
    '''
    Creates a plot of gene completeness. Plot is a bar chart with genes sorted from highest to lowest completeness on 
    the x-axis and percent complete on the y-axis. Completeness is calculated as percent of taxa_comp_df the gene is 
    present in.

    :param s: series containing completeness values for genes or taxa
    :type s: pd.Series
    :param plot_name: name of the plot file to be saved
    :type plot_name: str
    :param y_count: number of genes or taxa used to calculate completeness
    :type y_count: int
    :param genes: whether the completeness values are for genes or taxa, defaults to True
    :type genes: bool, optional
    '''

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


def get_md5(filename):
    '''
    Calculate the MD5 hash of a file.

    :param filename: path to the file to be hashed
    :type filename: str
    :return: MD5 hash of the file as a hexadecimal string
    :rtype: str
    '''
    # Open,close, read file and calculate MD5 on its contents
    with open(filename, 'rb') as file_to_check:
        # read contents of the file
        data = file_to_check.read()
        # pipe contents of the file through
        md5_returned = hashlib.md5(data).hexdigest()

    return md5_returned    


def backup(dfo):
    '''
    Creates a backup of the phylofisher database in the backups directory.

    :param dfo: path to the directory containing the phylofisher database
    :type dfo: str
    '''
    now = datetime.now().strftime('%d-%b-%Y_%H-%M-%S')

    if os.path.isdir(f'{dfo}/backups'):
        backups = os.listdir(f'{dfo}/backups')
        backups.sort(key=lambda date: datetime.strptime(date, '%d-%b-%Y_%H-%M-%S'), reverse=True)

        back_files = {}
        db_files = {}
        os.mkdir(f'{dfo}/backups/{now}')

        for root, dirs, files in os.walk(f'{dfo}/backups'):
            for file in files:
                back_files[file] = (os.path.join(root, file))

        db_files['phylofisher.db'] = f'{dfo}/phylofisher.db'

        for k, v in db_files.items():
            if k == 'phylofisher.db':
                dest = k
            else:
                dest = '/'.join(v.split('/')[-2:])

            
            shutil.copy(v, f'{dfo}/backups/{now}/{dest}')

    else:
        os.mkdir(f'{dfo}/backups')
        os.mkdir(f'{dfo}/backups/{now}')
        shutil.copy(f'{dfo}/phylofisher.db', f'{dfo}/backups/{now}/phylofisher.db')
