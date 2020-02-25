![PhyloFIsher](docs/_static/fisher.png)

## 1. Installation via Conda

1. Install conda via Miniconda or Anaconda:<br/>
https://docs.conda.io/projects/conda/en/latest/user-guide/install/

2. Prepare conda virtual environment:<br/>
 `conda create -n fisher`
 
3. Activate conda environment:<br/>
`conda activate fisher`

4. Add PhyloFisher Anaconda Cloud channel to your channels:<br/>
`conda config --append channels phylofisher`

5. Install phylofisher:<br/>
`conda install phylofisher`

- After you finish using PhyloFisher, use `conda deactivate` to deactivate the `fisher` conda virtual enviornment.

## 2. Custom dataset preparation
1. Create a directory that will contain the dataset.

2. Create a subdirectory of the dataset directory named `orthologs`
    - This directory contains FASTA formatted sequence files containing orthologous gene seqeunces.
        - A single organism should NOT be present more than once per gene. 
    - Name each file as such: _geneA_.fas.
    - Each file should contain one sequence per organism (named by some short name).
        - Example: >ParatetrGEN<br/>

3. If you have paralogs (if not continue to Step 4), create a subdirectry of the dataset 
    directory named `paralogs`: <br/>
    - This directory contains FASTA formatted sequence files containing paralogous gene seqeunces.
        - A single organism should be present no more than once per gene. 
    - Name each file as such: _geneA_\_paralogs.fas
    - FASTA sequence descriptions should be: organism short name, two dots, 'p' and 5 digits.
        - Example: >ParatetrGEN..p88706<br/>
    
4. Copy orthomcl folder from our dataset to your dataset directory.
    * `PUT COMMAND HERE`

##### Resulting Dataset Directory Structure
    dataset_dir/
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
` build_dataset.py [options]`
        
## 3. Usage

**Insert Work Flow Figure Here**

1. Make directory for your analysis and go there:<br/>
`mkdir <analysis_dir>`<br/>
`cd <analysis_dir>`<br/>

2. (Optional) Prepare config file in your folder:<br/>
`config.py [OPTIONS] -d <dataset_folder> -i <input_file.tsv>`<br/>
Use config.py --help for more options.<br/>

3. Run fisher.py for ortholog fishing:<br/>
`fisher.py [OPTIONS]`<br/>

4. Prepare *_stats folder for organism and gene selection:<br/>
`informant.py [OPTIONS] -i fasta/ --paralog_selection`<br/>

5. Apply your selection:<br/>
`fishing_net.py [OPTIONS] -i fasta/ -o <folder_with_selected_seqs>`<br/>

6. Trimming and single gene tree computations:<br/>
`trimming.py [OPTIONS] -d <dataset_folder> -i <folder_with_fasta_files> -s <suffix>`

7. Build phylogenomic matrix (Concatenation):<br/>
`forge.py [OPTIONS] --input <input_dir> `



## Utilities
#### examine_bipartitions.py<br/>
* Description
* Usage: `usage goes here`
    * Use `examine_bipartitions.py --help` to see all options
    
#### fast_site_removal.py<br/>
* Description
* Usage: `usage goes here`
    * Use `fast_site_removal.py --help` to see all options

#### aa_comp_calculator.py<br/>
* Description
* Usage: `usage goes here`
    * Use `aa_comp_calculator.py --help` to see all options

#### fast_tax_removal.py<br/>
* Description
* Usage: `usage goes here`
    * Use `fast_tax_removal.py --help` to see all options

#### genetic_code_check.py<br/>
* Description
* Usage: `usage goes here`
    * Use `genetic_code_check.py --help` to see all options

#### heteroevoloving_sites.py<br/>
* Description
* Usage: `usage goes here`
    * Use `heteroevoloving_sites.py --help` to see all options

#### SR4_class_recoder.py<br/>
* Description
* Usage: `usage goes here`
    * Use `SR4_class_recoder.py --help` to see all options

#### missing_data.py<br/>
* Description
* Usage: `usage goes here`
    * Use `missing_data.py --help` to see all options