import os
import subprocess
import pytest
import shutil

DATA_DIRECTORY = os.path.join(os.path.dirname(__file__), 'data')

def test_bipartition_examiner():
    '''
    Test bipartition_examiner.py
    '''
    bipar_data_dir = os.path.join(DATA_DIRECTORY, 'bipartition_examiner')
   
   # Test no db
    cmd = f'bipartition_examiner.py -b {bipar_data_dir}/bs_files.txt -g {bipar_data_dir}/no_db_groups.txt -o bipar_test_output --no_db'
    subprocess.run(cmd, shell=True)
    assert os.path.isfile('bipar_test_output/bipartition_examiner.tsv')
    assert os.path.isfile('bipar_test_output/bipartition_examiner.pdf')
    shutil.rmtree('bipar_test_output')
    
    # TODO: Add test for with db option