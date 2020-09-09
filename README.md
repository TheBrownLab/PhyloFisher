![PhyloFIsher](docs/TREE-PF-LOGO-WOtree.png)

Jump to:
   - [Installation](#Installation)
   - [Usage](#Usage)
   - [Custom Dataset Prep](#Custom-dataset-preparation)

<br/>

# Installation

### via conda (Recommended)
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

### via pip
1. Install PhyoFisher:<br/>
`pip install phylofisher`
2. Add the following line to your .bashrc or .bash_profile:<br/>
`export PATH="$HOME/bin:$PATH"`
3. Reload your .bashrc or .bash_profile<br/>
`source <.bashrc | .bash_profile>`<br/>
4. Install Non-Python Dependencies:<br/>
Linux: `install_deps.py`
macOS: `install_deps.py -gxx <path/to/g++>`

Notes:
- For OSX users, cd-hit requires gcc for compilation (to use Homebrew, see https://brewformulas.org/gcc).  

# Usage

**Insert Work Flow Figure Here**

### Main Workflow
1. **Prepare Input Metadata**
1. **config.py**
    * config.py allows the user to set to configure the analyis directory.
    * Usage: `config.py [OPTIONS] -d <database_folder> -i <input_metadata.tsv>`<br/>
1. **fisher.py**
    * Run fisher.py for ortholog fishing.
    - Usage: `fisher.py [OPTIONS]`<br/>
1. **informant.py**
    * Prepare *_stats folder for organism and gene selection.
    * Usage: `informant.py [OPTIONS] -i fasta/ --paralog_selection`<br/>
1. **fishing_net.py**
    * Apply your selection.
    * Usage: `fishing_net.py [OPTIONS] -i fasta/ -o <folder_with_selected_seqs>`<br/>
1. **sgt_constructor.py**
    * Length filtration, trimming and single gene tree computations.
    * Usage: `trimming.py [OPTIONS] -i <folder_with_fasta_files>`<br/>
1. **forest.py**
    * Tree visualization.
    * Usage: `forest.py [OPTIONS] --input <dir_containing_tress>`<br/>
1. **lumberjack.py**
    * Add new data to the starting dataset.
    * Usage: `INSERT USAGE`
1. **select_taxa.py** (OPTIONAL)
    * Select a subset of taxa from the database for final phylogenomic analyses.
    * Usage: `select_taxa.py [OPTIONS]`<br/>
1. **select_orthologs.py** (OPTIONAL)
    * Select a subset of orthologs from the database for final phylogenomic analyses.
    * Usage: `select_orthologs.py [OPTIONS]`<br/>
1. **prep_final_dataset.py**
    * Collect taxa and genes for final phylogenomic analyses.
    * Usage: `prep_final_dataset.py [OPTIONS]`<br/>
1. **forge.py**
    * Build phylogenomic matrix (Concatenation).
    * Usage: `forge.py [OPTIONS] --input <input_dir> `<br/>
 

***Use `<script>.py --help`  to see all options and their desctriptions.


### Utilities
* **aa_comp_calculator.py**
    * Calculation of amino acid composition and using euclidean distances to hierarchically 
    cluster these data, in order to examine if amino acid composition may contribute to the
     groupings that were inferred in phylogenies. This type of analysis was used in Brown 
     et al. 2018 as an example. 
    * Usage: `aa_comp_calculator.py [OPTIONS] -i <matrix>`<br/>
* **astral_runner.py**
    * Generates input files and infers a coalescent-based species tree 
    given a set of single ortholog trees and bootstrap trees using ASTRAL-III.
    * Usage: `a`<br/>
* **backup_restoration.py**
    * Restores the database from a backup.
    * Usage: `a`<br/>
* **bipartition_examiner.py**
    * Calculates the observed occurrences of clades of interest in bootstrap trees.
    * Usage: `bipartition_examiner.py [OPTIONS] -i /path/to/input/`<br/>
* **fast_site_removal.py**
    * The fastest evolving sites are expected to be the most prone to phylogenetic 
    signal saturation and systematic model misspecification in phylogenomic analyses. 
    This tool will remove the fastest evolving sites within the phylogenomic supermatrix 
    in a stepwise fashion, leading to a user defined set of new matrices with these 
    sites removed.
    * Usage: `fast_site_removal.py [OPTIONS] -i /path/to/input/`<br/>
* **fast_tax_removal.py**
    * Removes the fastest evolving taxa, based on branch length. This tool will remove the
     fastest evolving taxa within the phylogenomic supermatrix in a stepwise fashion, 
     leading to a user defined set of new matrices with these taxa removed.
    * Usage: `fast_tax_removal.py [OPTIONS] -t <tree> -m <matrix> -i <num_of_iterations>`<br/>
* **genetic_code_check.py**
    * Script for fast genetic code analysis.
    * Usage: `genetic_code.py [OPTIONS] -i, --input file.fas -q, --queries Allomacr,Mantplas,...`<br/>
* **heterotachy.py**
    * Within-site rate variation (heterotachy) (Lopez Casane & Philippe, 2002) has been
     shown to cause artificial relationships in molecular phylogenetic reconstruction 
     (Inagaki et al., 2004).  This tool will remove the most heterotacheous sites within 
     a phylogenomic supermatrix in a stepwise fashion, leading to a user defined set of 
     new matrices with these sites removed.
    * Usage: `heterotachy.py -t path/to/tree -m path/to/matrix [OPTIONS]`<br/>
* **mammal_modeler.py**
    * Generates a MAMMaL site heterogeneous model from a user input tree and supermatrix 
    with estimated frequencies for a user defined number of classes using the methods 
    described in ADD CITATION.
    * Usage: `h`<br/>
* **purge.py**
    * This tool is used for deleting taxa and/or taxonomic groups from the dataset permanently. 
    * Usage: `h`<br/>
* **random_sample_iteration.py**
    * This tool randomly resamples the gene set into a set of new matrices that are subsamples 
    of the super matrix. It constructs supermatrices from randomly sampled genes with user 
    defined options such as the confidence interval sampling all genes in a random fashion and 
    the percentage of subsampling a user requires per sampling step. This method was used in Brown 
    et al. 2018 as an example. 
    * Usage: `h`<br/>
* **rtc_binner.py**
    * Calculates the relative tree certainty score (RTC) in RAxML of each single ortholog 
    tree and bins them based on their RTC scoring into top 25%, 50%, and top 75% sets. 
    Supermatrices are constructed from these bins of orthologs. 
    * Usage: `S`<br/>
* **SR4_class_recoder.py**
    * To minimize phylogenetic saturation this tool recodes input supermatrix into the 
    four-character state scheme of SR4 (Susko and Roger 2007), based on  amino acid classification.
    * Usage: `SR4_class_recoder.py [OPTIONS] -i <matrix>`<br/>
* **taxon_collapser.py**
    * Allows users to combine multiple operational taxonomic units into one single taxon. 
    For example if a user has multiple single cell libraries from a taxon or multiple strains 
    of the same species (or genus etc.), a user may decide to  collapse all these 
    strains/libraries into a single taxon. 
    * Usage: `S`<br/>

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
        
