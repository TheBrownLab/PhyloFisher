![PhyloFIsher](docs/_static/fisher.png)

## 1. Installation
### 1.1 via Conda

- Install conda via Miniconda or Anaconda:<br/>
https://docs.conda.io/projects/conda/en/latest/user-guide/install/

- Prepare conda virtual environment:<br/>
 `conda create -n fisher`
 
- Activate conda environment:<br/>
`conda activate fisher`

- Install phylofisher:<br/>
`conda install -c _channel_ phylofisher`

- After you finish using PhyloFisher use `conda deactivate` to deactivate `fisher`

## 2. Custom dataset preparation
- Create a directory with subdirectory _orthologs_ (_genename_.fas) with one sequence per organism 
(named by some short name)

- If you already have paralogs, they should be named like this:
    - ParatetrGEN..p88706 
        - Organism short name, two dots, 'p' and 5 digits
    - Stored in 'paralogs' folder like _genename_\_paralogs.fas
    
- Copy orthomcl folder from our dataset to your folder

##### Dataset Directory Structure
    dataset_dir/
        orthologs/
            geneA.fas
            geneB.fas
            ...
            geneZ.fas
        paralogs/
            geneA_paralogs.fas
            geneB_paralogs.fas
            ...
            geneZ_paralogs.fas
        orthomcl/

- Run _build_dataset.py_<br/>
` build_dataset.py [options]`
        
## 3. Usage

- Make directory for your analysis and go there:<br/>
`mkdir <analysis_dir>`<br/>
`cd <analysis_dir>`<br/>

- Prepare config file in your folder:<br/>
`config.py -d <dataset_folder> -i <input_file.tsv> [OPTIONS]`<br/>
Use config.py --help for more options.<br/>

- Run fisher.py for ortholog fishing:<br/>
`fisher.py [OPTIONS]`<br/>

- Prepare *_stats folder for organism and gene selection:<br/>
`informant.py -i fasta/ --paralog_selection`<br/>

- Apply your selection:<br/>
`fishing_net.py -i fasta/ -o <folder_with_selected_seqs>`<br/>

- Trimming and trees computing<br/>
`trimming.py -d <dataset_folder> -i <folder_with_fasta_files> -s <sufix>`



## Utilities

- add_aln_length.py
    - description
    - 
- AgentEllie.py
- AgentEllie2.py
- CalculateAAComposition.py
- corrected_translation.py
- fast_tax_removal.py
- genetic_code_check.py
- heteroevoloving_sites.py
- len_filter.py
- len_filter2.py
- no_gap_stops.py
- pre_trimal.py
- Recode_Phylip_SR4classes.py