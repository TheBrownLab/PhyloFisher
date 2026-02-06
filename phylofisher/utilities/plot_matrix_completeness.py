#!/usr/bin/env python
import os
import textwrap
from Bio import SeqIO
import matplotlib.pyplot as plt
from phylofisher import help_formatter

def get_taxa_from_list(taxa_file):
    '''
    Extracts taxa from a given text file, maintaining the order from the file.

    :param taxa_file: Path to the text file containing taxa names (one per line).
    :type taxa_file: str
    :return: List of taxa names from the file in order.
    :rtype: list
    '''
    taxa_list = []
    with open(taxa_file, 'r') as f:
        for line in f:
            taxon = line.strip()
            if taxon:  # Skip empty lines
                taxa_list.append(taxon)

    # Reverse the order for plotting (to maintain visual consistency with tree plotting)
    taxa_list.reverse()

    return taxa_list


def get_gene_completeness(occupancy_file, taxa):
    '''
    Calculates the gene-wise completeness for the given taxa from an occupancy TSV file.

    :param occupancy_file: Path to the occupancy TSV file.
    :type occupancy_file: str
    :param taxa: List of taxa names to check completeness for.
    :type taxa: list
    :return: Dictionary with taxa as keys and their gene completeness as values.
    :rtype: dict
    '''
    gene_completeness = {}
    
    with open(occupancy_file, 'r') as f:
        # Read header line to get gene names
        header = f.readline().strip().split('\t')
        genes = header[1:]  # Skip the first empty column
        
        # Read data lines
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue
            
            taxon = parts[0]
            if taxon in taxa:
                # Count presence/absence (1s and 0s)
                gene_presence = [int(x) for x in parts[1:]]
                completeness = sum(gene_presence) / len(gene_presence) if gene_presence else 0
                gene_completeness[taxon] = completeness
    
    # Check if all taxa were found
    for taxon in taxa:
        if taxon not in gene_completeness:
            raise ValueError(f'Taxon {taxon} not found in the provided occupancy file.')
    
    return gene_completeness


def get_completeness(matrix_file, taxa):
    '''
    Calculates the site-wise completeness of the matrix for the given taxa.

    :param matrix_file: Path to the matrix file.
    :type matrix_file: str
    :param taxa: List of taxa names to check completeness for.
    :type taxa: list
    :return: Dictionary with taxa as keys and their site-wise completeness as values.
    :rtype: dict
    '''
    completeness = {}
    with open(matrix_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id in taxa:
                completeness[record.id] = 1 - (record.seq.count('-') / len(record.seq))
            else:
                raise ValueError(f'Taxon {record.id} not found in the provided taxa list.')
    return completeness


def plot_completeness(completeness, gene_completeness, taxa_order, output_filename='completeness.pdf'):
    '''
    Plots the site-wise and gene-wise completeness of the taxa in the specified order.

    :param completeness: Dictionary with taxa as keys and their site-wise completeness as values.
    :type completeness: dict
    :param gene_completeness: Dictionary with taxa as keys and their gene-wise completeness as values.
    :type gene_completeness: dict
    :param taxa_order: List of taxa names in the desired plot order.
    :type taxa_order: list
    :param output_filename: Name of the output PDF file.
    :type output_filename: str
    '''

    # Use the provided taxa order instead of dictionary keys
    taxa = taxa_order
    seq_values = [completeness[taxon] for taxon in taxa_order]
    gene_values = [gene_completeness[taxon] for taxon in taxa_order]

    # Create positions for the bars
    y_pos = range(len(taxa))
    bar_height = 0.35

    plt.figure(figsize=(4.5, len(taxa) * 0.25))
    
    # Add 0%, 50%, and 100% lines behind the bars
    plt.axvline(x=0, color='gray', linestyle='-', alpha=0.7, zorder=1, linewidth=0.3)
    plt.axvline(x=0.5, color='gray', linestyle='-', alpha=0.7, zorder=1, linewidth=0.3)
    plt.axvline(x=1, color='gray', linestyle='-', alpha=0.7, zorder=1, linewidth=0.3)

    # Create horizontal bars
    gene_bars = plt.barh([y + bar_height/2 for y in y_pos], gene_values, 
                         bar_height, color='#74a446', zorder=2, label='Gene-wise')
    seq_bars = plt.barh([y - bar_height/2 for y in y_pos], seq_values, 
                        bar_height, color='#446074', zorder=2, label='Site-wise')
    
    plt.xlim(0, 1)
    plt.ylim(-0.5, len(taxa) - 0.5)
    plt.xticks([0, 0.5, 1], ['0%', '50%', '100%'], rotation=270)
    plt.yticks(y_pos, taxa)
    
    # Add legend
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Remove border and customize grid
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.tick_params(axis='x', length=0)
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight', pad_inches=0.1)


if __name__ == '__main__':
    description = 'Plots site-wise and gene-wise completeness for taxa in a given taxa list.'
    parser, optional, required = help_formatter.initialize_argparse(name='plot_matrix_completeness.py',
                                                                    desc=description,
                                                                    usage='plot_matrix_completeness.py '
                                                                          '[OPTIONS] -i <matrix> -l <taxa_list> -g <gene_occupancy>')

    required.add_argument('-i', '--input', required=True, type=str, metavar='matrix',
                          help=textwrap.dedent("""\
                          Path to input matrix for analysis.
                          """))
    required.add_argument('-l', '--taxa_list', required=True, type=str, metavar='taxa_list',
                          help=textwrap.dedent("""\
                          Path to text file containing taxa names (one per line).
                          """))
    required.add_argument('-g', '--gene_occupancy', required=True, type=str, metavar='gene_occupancy',
                          help=textwrap.dedent("""\
                          Path to gene occupancy TSV file.
                          """))

    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=False)

    taxa_file = os.path.abspath(args.taxa_list)
    matrix_file = os.path.abspath(args.input)
    occupancy_file = os.path.abspath(args.gene_occupancy)
    basename = '.'.join(os.path.basename(matrix_file).split('.')[:-1])
    output_filename = f'{basename}_completeness.pdf'
    
    taxa = get_taxa_from_list(taxa_file)
    completeness = get_completeness(matrix_file, taxa)
    gene_completeness = get_gene_completeness(occupancy_file, taxa)

    plot_completeness(completeness, gene_completeness, taxa, output_filename)