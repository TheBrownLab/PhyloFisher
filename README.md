![PhyloFIsher](docs/_static/fisher.png)

Jump to:
   - [Installation](#Installation-via-Conda)
   - [Usage](#Usage)
   - [Custom Dataset Prep](#Custom-dataset-preparation)

<br/>

# Installation

###via conda (Recommended)
1. Install Anaconda:<br/>
https://docs.conda.io/projects/conda/en/latest/user-guide/install/
2. Prepare a conda virtual environment:<br/>
 `conda create -n fisher`
3. Activate the conda environment:<br/>
`conda activate fisher`
4. Add the Bioconda, Conda-Forge, & PhyloFisher Anaconda Cloud channels to your channels:<br/>
`conda config --append channels phylofisher`<br/>
`conda config --append channels bioconda`<br/>
`conda config --append channels conda-forge`<br/>
5. Install PhyloFisher:<br/>
`conda install phylofisher`

Notes:
- After you finish using PhyloFisher, use `conda deactivate` to deactivate the `fisher` conda virtual enviornment.
<br/>

###via pip
1. Install PhyoFisher:<br/>
`pip install phylofisher`
2. Install Non-Python Dependencies:<br/>
`install_deps.py [OPTIONS]`

Notes:
- For OSX users, cd-hit requires gcc for compilation (to use Homebrew, see https://brewformulas.org/gcc).  

# Usage

**Insert Work Flow Figure Here**

### Main Workflow
1. **config.py**
    - config.py allows the user to set to configure the analyis directory <br/>
    - Usage: `config.py [OPTIONS] -d <dataset_folder> -i <input_metadata.tsv>`<br/>
2. **fisher.py**
    - Run fisher.py for ortholog fishing:<br/>
    - Usage: `fisher.py [OPTIONS]`<br/>
3. **informant.py**
    - Prepare *_stats folder for organism and gene selection:<br/>
    - Usage: `informant.py [OPTIONS] -i fasta/ --paralog_selection`<br/>
4. **fishing_net.py**
    - Apply your selection:<br/>
    - Usage: `fishing_net.py [OPTIONS] -i fasta/ -o <folder_with_selected_seqs>`<br/>
5. **single_gene_tree_constructor.py**
    - Trimming and single gene tree computations:<br/>
    - Usage: `trimming.py [OPTIONS] -i <folder_with_fasta_files>`<br/>
6. **forest.py**
    - Tree visualization:<br/>
    - Usage: `forest.py [OPTIONS] --input <dir_containing_tress>`<br/>
7. **forge.py**
    - Build phylogenomic matrix (Concatenation):<br/>
    - Usage: `forge.py [OPTIONS] --input <input_dir> `<br/>
8. **lumberjack.py**
    - Add new data to the dataset:<br/>
    - Usage: `lumberjack.py [OPTIONS] --input <path/to/tsvs>` <br/>

***Use `<script>.py --help`  to see all options and their desctriptions.


### Utilities
* **bipartition_examiner.py**
    * Description
    * Usage: `bipartition_examiner.py [OPTIONS] -i /path/to/input/`<br/>
* **fast_site_removal.py**
    * Description
    * Usage: `usage goes here`<br/>
* **aa_comp_calculator.py**
    * Calculates amino acid compostition of the supplied super matrix
    * Usage: `aa_comp_calculator.py [OPTIONS] -i <matrix>`<br/>
* **fast_tax_removal.py**
    * Description
    * Usage: `usage goes here`<br/>
* **genetic_code_check.py**
    * Script for fast genetic code analysis.
    * Usage: `genetic_code.py [OPTIONS] -i, --input file.fas -q, --queries Allomacr,Mantplas,...`<br/>
* **heteroevoloving_sites.py**
    * Description
    * Usage: `usage goes here`<br/>
* **SR4_class_recoder.py**
    * Recodes input matrix based on SR4 amino acid classifications.
    * Usage: `SR4_class_recoder.py [OPTIONS] -i <matrix>`<br/>
* **missing_data.py**
    * Subsets gene based on gene completeness.
    * Usage: `missing_data.py [OPTIONS] -i <input_dir> -m <metadata> {-n <gene_number> |
     -c <percent_complete>}`<br/>

***Use `<script>.py --help`  to see all options and their desctriptions.
 
 
      
# Custom dataset preparation
1. Create a directory that will contain the dataset and go there:<br/>
`mkdir dataset`<br/>
`cd dataset`<br/><br/>
2. Create a subdirectory of the dataset directory named orthologs (`mkdir orthologs`):
    - This directory contains FASTA formatted sequence files containing orthologous gene seqeunces.
        - A single organism should NOT be present more than once per gene. 
    - Name each file as such: _geneA_.fas.
    - Each file should contain one sequence per organism (named by some short name).
        - Example: >ParatetrGEN<br/><br/>
3. If you have paralogs (if not continue to Step 4), create a subdirectry of the dataset 
    directory named paralogs (`mkdir paralogs`):
    - This directory contains FASTA formatted sequence files containing paralogous gene seqeunces.
    - Name each file as such: _geneA_\_paralogs.fas
    - FASTA sequence descriptions should be: organism short name, two dots, 'p' and 5 digits.
        - Example: >ParatetrGEN..p88706<br/><br/>
4. Copy the OrthoMCL dataset to your dataset directory.
    * `PUT COMMAND HERE`
    
##### Resulting Dataset Directory Structure
    dataset/
        orthologs/
            geneA.fas
            geneB.fas
            ...
            geneZ.fas
        paralogs***/
            geneA_paralogs.fas
            geneB_paralogs.fas
            ...
            geneZ_paralogs.fas
        orthomcl/
*** only present if you already have paralogs

* Run _build_dataset.py_<br/>
` build_dataset.py [OPTIONS]`
    * Run `build_dataset.py --help` to see all options and their descriptions.
        