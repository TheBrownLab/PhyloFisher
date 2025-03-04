import os
import sqlite3
import pandas as pd
import shutil
from Bio import SeqIO
from phylofisher import help_formatter

def connect_to_db(db_file):
    '''
    Establish a connection to the SQLite database.

    :param db_file: path to the SQLite database file
    :type db_file: str
    :raises Exception: sqlite3.Error
    :return: sqlite3 connection object
    :rtype: sqlite3.Connection
    '''
    try:
        return sqlite3.connect(db_file)
    except sqlite3.Error as e:
        raise Exception(f'Database connection error: {e}')

def create_columns(cursor, genes):
    '''
    Create tables in the SQLite database for taxonomies, metadata, and sequences.

    :param cursor: SQLite database cursor
    :type cursor: sqlite3.Cursor
    :param genes: list of genes
    :type genes: list
    '''
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS taxonomies (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            taxonomy VARCHAR(255) UNIQUE,
            color TEXT
        );
    ''')

    cursor.execute('''
        CREATE INDEX IF NOT EXISTS idx_taxonomy ON taxonomies (taxonomy);
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS metadata (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            short_name VARCHAR(255) UNIQUE,
            long_name TEXT,
            higher_taxonomy_id INTEGER,
            lower_taxonomy_id INTEGER,
            data_type TEXT NOT NULL CHECK (data_type IN ('Genomic', 'Transcriptomic', 'EST', 'Collapsed')),
            source TEXT,
            FOREIGN KEY (higher_taxonomy_id) REFERENCES taxonomies(id),
            FOREIGN KEY (lower_taxonomy_id) REFERENCES taxonomies(id)
        );
    ''')

    cursor.execute('''
        CREATE INDEX IF NOT EXISTS idx_short_name ON metadata (short_name);
    ''')
    
    for gene in genes:
        table_name = f'"{gene}"'

        cursor.execute(f'''
            CREATE TABLE IF NOT EXISTS {table_name} (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                header TEXT,
                sequence TEXT,
                is_paralog BOOLEAN,
                metadata_id INTEGER,
                FOREIGN KEY (metadata_id) REFERENCES metadata(id)
            );
        ''')

        cursor.execute(f'''
            CREATE INDEX IF NOT EXISTS idx_header ON {table_name} (header);
        ''')

def load_taxonomies(cursor, metadata_tsv_file, tree_colors_file):
    '''
    Parse the metadata.tsv and tree_colors.tsv files to get the taxonomies and their associated colors,
    and load them into the taxonomies table in the database.

    :param cursor: SQLite database cursor
    :type cursor: sqlite3.Cursor
    :param metadata_tsv_file: path to the metadata.tsv file
    :type metadata_tsv_file: str
    :param tree_colors_file: path to the tree_colors.tsv file
    :type tree_colors_file: str
    '''
    # Parse the metadata and tree_colors TSV files using pandas
    try:
        with open(metadata_tsv_file, 'r') as metadata_tsv, open(tree_colors_file, 'r') as tree_colors:
            metadata_df = pd.read_csv(metadata_tsv, sep='\t')
            tree_colors_df = pd.read_csv(tree_colors, sep='\t')
            
            # Extract all unique taxonomies from the metadata file
            taxonomies = set(metadata_df['Higher Taxonomy'].to_list()).union(set(metadata_df['Lower Taxonomy'].to_list()))

            # Create a dictionary to map taxonomies to their respective colors
            taxonomy_colors = {}
            for taxonomy in taxonomies:
                # Get the color from the tree_colors file, or assign None if not found
                color = tree_colors_df.loc[tree_colors_df['Taxonomy'] == taxonomy, 'Color'].values
                taxonomy_colors[taxonomy] = color[0] if color.size > 0 else None

            # Insert taxonomies and their colors into the database
            for taxonomy, color in taxonomy_colors.items():
                cursor.execute('''
                    INSERT INTO taxonomies (taxonomy, color)
                    VALUES (?, ?)
                    ON CONFLICT(taxonomy) DO UPDATE SET color=excluded.color;
                ''', (taxonomy, color))

    except Exception as e:
        raise Exception(f"Error while loading taxonomies: {e}")

def load_metadata(cursor, metadata_tsv_file):
    '''
    Load metadata from a TSV file into the metadata table using pandas.

    :param cursor: SQLite database cursor
    :type cursor: sqlite3.Cursor
    :param metadata_tsv_file: Path to the metadata TSV file
    :type metadata_tsv_file: str
    '''
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(metadata_tsv_file, sep='\t')

    # Create a list to hold the values to be inserted
    insert_values = []

    for _, row in df.iterrows():
        # Retrieve the higher and lower taxonomy IDs
        higher_taxonomy = row['Higher Taxonomy']
        lower_taxonomy = row['Lower Taxonomy']

        # Query the taxonomies table to get the IDs
        cursor.execute('SELECT id FROM taxonomies WHERE taxonomy = ?', (higher_taxonomy,))
        higher_taxonomy_id = cursor.fetchone()
        higher_taxonomy_id = higher_taxonomy_id[0] if higher_taxonomy_id else None

        cursor.execute('SELECT id FROM taxonomies WHERE taxonomy = ?', (lower_taxonomy,))
        lower_taxonomy_id = cursor.fetchone()
        lower_taxonomy_id = lower_taxonomy_id[0] if lower_taxonomy_id else None

        # Prepare the data for insertion
        insert_values.append((
            row['Unique ID'],  # short_name
            row['Long Name'],   # long_name
            higher_taxonomy_id, # higher_taxonomy_id
            lower_taxonomy_id,  # lower_taxonomy_id
            row['Data Type'],   # data_type
            row['Source']       # source
        ))

    # Insert data into metadata table
    cursor.executemany('''
        INSERT INTO metadata (short_name, long_name, higher_taxonomy_id, lower_taxonomy_id, data_type, source)
        VALUES (?, ?, ?, ?, ?, ?)
    ''', insert_values)

def load_sequences(cursor, metadata_tsv_file, old_database_dir, genes):
    '''
    Load sequences from the old database directory into the sequences table in the new database.

    :param cursor: SQLite database cursor
    :type cursor: sqlite3.Cursor
    :param metadata_tsv_file: path to the metadata.tsv file
    :type metadata_tsv_file: str
    :param old_database_dir: path to the old database directory
    :type old_database_dir: str
    :param genes: list of genes
    :type genes: list
    '''
    for gene in genes:
        table_name = f'"{gene}"'  
        insert_values = []  

        # Get path to ortholog and paralog files
        paralog_file = os.path.join(old_database_dir, 'paralogs', f'{gene}_paralogs.fas')
        ortholog_file = os.path.join(old_database_dir, 'orthologs', f'{gene}.fas')

        # Check if files exist before attempting to parse
        if not os.path.exists(ortholog_file):
            print(f"Warning: Missing ortholog file {ortholog_file}")
        else:
            for record in SeqIO.parse(ortholog_file, 'fasta'):
                cursor.execute('SELECT id FROM metadata WHERE short_name = ?', (record.id,))
                metadata_id = cursor.fetchone()
                if metadata_id is None:
                    print(f"Warning: No metadata entry found for {record.id}")
                    continue
                insert_values.append((record.id, str(record.seq), False, metadata_id[0]))

        # Check if files exist before attempting to parse
        if not os.path.exists(paralog_file):
            print(f"Warning: Missing paralog file {paralog_file}")
        else:
            for record in SeqIO.parse(paralog_file, 'fasta'):
                short_name = record.id.split('..')[0]
                cursor.execute('SELECT id FROM metadata WHERE short_name = ?', (short_name,))
                metadata_id = cursor.fetchone()
                if metadata_id is None:
                    print(f"Warning: No metadata entry found for {short_name}")
                    continue
                insert_values.append((short_name, str(record.seq), True, metadata_id[0]))

        # Execute insertion only if data exists
        if insert_values:
            cursor.executemany(f'''
                INSERT INTO {table_name} (header, sequence, is_paralog, metadata_id)
                VALUES (?, ?, ?, ?)
            ''', insert_values)


def cp_everything_else(old_database_dir, output_dir):
    '''
    Copy orthomcl, proteomes, datasetdb, and profiles directories to the new database directory.

    :param old_database_dir: path to the old database directory
    :type old_database_dir: str
    :param output_dir: path to the new database directory
    :type output_dir: str
    '''
    required_dirs = ['orthomcl', 'proteomes', 'datasetdb', 'profiles']
    for sub_dir in required_dirs:
        src = os.path.join(old_database_dir, sub_dir)
        dest = os.path.join(output_dir, sub_dir)
        if os.path.exists(src):
            try:
                shutil.copytree(src, dest, dirs_exist_ok=True)
            except Exception as e:
                print(f'Error copying {sub_dir}: {e}')


if __name__ == '__main__':
    description = 'Script for converting an existing PhyloFisher dataset into a SQLite database.'
    parser, optional, required = help_formatter.initialize_argparse(name='dataset_to_database.py', desc=description, usage='dataset_to_database.py [OPTIONS]')
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=True, out_dir=True)
    
    # Set paths
    old_database_dir = args.input
    metadata_tsv_file = os.path.join(old_database_dir, 'metadata.tsv')
    tree_colors_tsv_file = os.path.join(old_database_dir, 'tree_colors.tsv')
    database_file = os.path.join(args.output, 'phylofisher.db')

    # Get genes from the orthologs directory. Important for creating the database schema for custom databases.
    genes = [f.split('.')[0] for f in os.listdir(os.path.join(old_database_dir, 'orthologs')) if f.endswith('.fas')]

    # Create the new database directory
    os.makedirs(args.output, exist_ok=True)

    # Initialize database connection
    conn = connect_to_db(database_file)
    cursor = conn.cursor()
    create_columns(cursor, genes)

    # Insert taxonomy data into the database
    load_taxonomies(cursor, metadata_tsv_file, tree_colors_tsv_file)

    # Insert metadata into the database
    load_metadata(cursor, metadata_tsv_file)

    # Insert sequences into the database
    load_sequences(cursor, metadata_tsv_file, old_database_dir, genes)

    # Commit and close database connection
    conn.commit()
    conn.close()

    # Copy other dataset components
    cp_everything_else(old_database_dir, args.output)
