import os
import sqlite3
import csv
import subprocess
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

def create_table(cursor, table_name, columns):
    '''
    Create a table if it doesn't exist.

    :param cursor: connection cursor
    :type cursor: str
    :param table_name: table name
    :type table_name: str
    :param columns: columns for table
    :type columns: str
    '''
    columns_def = ', '.join([f'`{col}` TEXT' for col in columns])
    query = f'''
    CREATE TABLE IF NOT EXISTS `{table_name}` (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        {columns_def}
    );
    '''
    cursor.execute(query)

def load_sequences_into_db(ortholog_dir, paralog_dir, db_file):
    '''
    Load FASTA sequences into the database.

    :param ortholog_dir: path to the ortholog FASTA files
    :type ortholog_dir: str
    :param paralog_dir: path to the paralog FASTA files
    :type paralog_dir: str
    :param db_file: path to the SQLite database file
    :type db_file: str
    '''
    conn = None
    try:
        conn = connect_to_db(db_file)
        cursor = conn.cursor()

        # Combine ortholog and paralog directories
        fasta_files = [
            (os.path.join(paralog_dir, f), 'paralog') for f in os.listdir(paralog_dir)
        ] + [
            (os.path.join(ortholog_dir, f), 'ortholog') for f in os.listdir(ortholog_dir)
        ]

        for fasta_file, ortholog_paralog in fasta_files:
            table_name = os.path.splitext(os.path.basename(fasta_file))[0].split('_')[0]
            create_table(cursor, table_name, ['header', 'sequence', 'ortholog_paralog'])

            # Parse and insert records
            data = []
            for record in SeqIO.parse(fasta_file, 'fasta'):
                record_id = record.id.split('..p')[0] if '..p' in record.id else record.id
                data.append((record_id, str(record.seq), ortholog_paralog))

            placeholders = ', '.join(['?' for _ in range(3)])
            insert_query = f'INSERT INTO `{table_name}` (header, sequence, ortholog_paralog) VALUES ({placeholders})'
            cursor.executemany(insert_query, data)

        conn.commit()
    except Exception as e:
        print(f'Error while loading sequences: {e}')
    finally:
        if conn:
            conn.close()

def load_tsv_into_db(tsv_file, db_file, table_name):
    '''
    Load TSV data into the database.

    :param tsv_file: path to the TSV file
    :type tsv_file: str
    :param db_file: path to the SQLite database file
    :type db_file: str
    :param table_name: SQLite table name
    :type table_name: str
    '''
    conn = None
    try:
        conn = connect_to_db(db_file)
        cursor = conn.cursor()

        # Read TSV data
        with open(tsv_file, 'r', newline='', encoding='utf-8') as file:
            reader = csv.reader(file, delimiter='\t')
            # Extract headers
            headers = next(reader)  

            create_table(cursor, table_name, headers)
            data = [row for row in reader]

            placeholders = ', '.join(['?' for _ in headers])
            insert_query = f'INSERT INTO `{table_name}` ({", ".join(f"`{col}`" for col in headers)}) VALUES ({placeholders})'
            cursor.executemany(insert_query, data)

        conn.commit()
    except Exception as e:
        print(f'Error while loading TSV: {e}')
    finally:
        if conn:
            conn.close()

def cp_everything_else(old_database_dir):
    '''
    Copies orthomcl, proteomes, datasetdb, and profiles directories to the new database directory.

    :param old_database_dir: path to the old database directory
    :type old_database_dir: str
    '''
    subprocess.run(f'cp -r {old_database_dir}/orthomcl {args.output}', shell=True, check=True)
    subprocess.run(f'cp -r {old_database_dir}/proteomes {args.output}', shell=True, check=True)
    subprocess.run(f'cp -r {old_database_dir}/datasetdb {args.output}', shell=True, check=True)
    subprocess.run(f'cp -r {old_database_dir}/profiles {args.output}', shell=True, check=True)


if __name__ == '__main__':
    description = 'Script for coverting a existing PhyloFisher dataset into a SQLite database.'
    parser, optional, required = help_formatter.initialize_argparse(name='dataset_to_database.py',
                                                                    desc=description,
                                                                    usage='dataset_to_database.py [OPTIONS]')
    args = help_formatter.get_args(parser, optional, required, pre_suf=False, inp_dir=True, out_dir=True)

    # Set paths
    old_database_dir = args.input
    ortholog_dir = os.path.join(old_database_dir, 'orthologs')
    paralog_dir = os.path.join(old_database_dir, 'paralogs')
    metadata_tsv_file = os.path.join(old_database_dir, 'metadata.tsv')
    tree_colors_tsv_file = os.path.join(old_database_dir, 'tree_colors.tsv')
    database_file = os.path.join(args.output, 'phylofisher.db')

    # Makes new database directory
    os.mkdir(args.output)

    # Load data into database
    load_sequences_into_db(ortholog_dir, paralog_dir, database_file)
    load_tsv_into_db(metadata_tsv_file, database_file, 'metadata')
    load_tsv_into_db(tree_colors_tsv_file, database_file, 'tree_colors')
    cp_everything_else(old_database_dir)