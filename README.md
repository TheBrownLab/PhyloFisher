![PhyloFIsher](fisher.png)

## 1. Installation
### 1.1 Installation with conda

- Install conda via Miniconda or Anaconda:<br/>
https://docs.conda.io/projects/conda/en/latest/user-guide/install/
- Prepare conda virtual environment:<br/>
 `conda create -n <your_env_name>`
- Activate phylo_fisher environment:<br/>
`conda activate <your_env_name>`
- Install phylofisher:<br/>
'conda install -c test_ phylofisher`
- After you finish using PhyloFisher use `conda deactivate` to deactivate `<your_env_name>`

## 2. Custom dataset preparation
- create a directory with subdirectory _orthologs_ (gene_name.fas) with one sequence per organism (named by some short name)
- ff you already have paralogs, they should be named like this:
ParatetrGEN..p88706 (short name or your organism, two dots, 'p' and 5 digits) and stored in 'paralogs' folder like genename_paralogs.fas
- copy orthomcl folder from our dataset to your folder
- run _build_dataset.py_
` build_dataset.py [options]`


## 3. Usage

- make directory for your analysis and go there:<br/>
`mkdir <analysis_dir>`<br/>
`cd <analysis_dir>`<br/>
- prepare config file in your folder:<br/>
`config.py -d <dataset_folder> -i <input_file.tsv> [OPTIONS]`<br/>
Use config.py --help for more options.<br/>
- run fisher.py for ortholog fishing:<br/>
`fisher.py [OPTIONS]`<br/>
- prepare *_stats folder for organism and gene selection:<br/>
`informant.py -i fasta/ --paralog_selection`<br/>
- apply your selection:<br/>
`fishing_net.py -i fasta/ -o <folder_with_selected_seqs>`<br/>
- trimming and trees computing<br/>
`trimming.py -d <dataset_folder> -i <folder_with_fasta_files> -s <sufix>`
