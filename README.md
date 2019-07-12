#Phylofisher
        
## 1. Installation
### 1.1 Installation with conda

- install conda via Miniconda or Anaconda:<br/>
https://docs.conda.io/projects/conda/en/latest/user-guide/install/
- prepare conda virtual environment:<br/>
 `conda env create -f fisher_env.yml`
- activate phylo_fisher environment:<br/>
`conda activate phylo_fisher`
- intall phylofisher:<br/>
`pip install .`
- after you finish using PhyloFisher use `conda deactivate` to deactivate phylo_fisher env

## 2. Usage

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
`NOTHING HERE SO FAR`
