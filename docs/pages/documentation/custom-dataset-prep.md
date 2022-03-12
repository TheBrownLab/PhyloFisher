---
layout: default
title: Custom Dataset Preparation
parent: Documentation
nav_order: 4
permalink: /documentation/custom-dataset-prep
---
# Custom Dataset Preparation

1. **Retrieve the recommended PhyloFisher directory structure, OrthoMCL v. 5.0 database, example
metadata.tsv, and tree_colors.tsv file via wget** <br/>
`wget https://ndownloader.figshare.com/files/29093325`
2. **Uncompress the file 29093325** <br/>
`tar -xzvf 29093325`
3. **Move into the directory PhyloFisher_FOR_CUSTOM_DATABASE** <br/>
`cd PhyloFisher_FOR_CUSTOM_DATABASE`
4. **Place your single gene files of orthologs in the directory /database/orthologs/**
    * The ortholog files must be in fasta format.
    * Each ortholog file must be named with the following convention {gene_name}.fas.
        * Ex. RPL7.fas
    * Each individual taxon should have a Unique ID as the header in all ortholog files. This Unique
        ID must be the the same in all ortholog files.
    * Each taxon can be present only once in each ortholog file.

5. **Place files of known paralogs for each gene in the directory /database/paralogs/ (OPTIONAL)**
    * Each gene file must be named with the following convention {gene_name}_paralogs.fas.
        * Ex. RPL7_paralogs.fas
    * Each individual taxon should have a Unique ID as the header in all paralog files. This Unique ID
        must be the the same in all paralog files and the corresponding ortholog files.
    * Each taxon can be present more than once in each paralog file.
    
6. **Place the complete proteome of each taxon present in the ortholog files in /database/proteomes**
    * All proteomes must be in fasta format
    * All proteomes must be tar and gzipped and follow the naming convention {Unique ID}.tar.gz

7. **Fill out the metadata.tsv file found in PhyloFisher_FOR_CUSTOM_DATABASE/database** </br>
   Detailed instructions on preparing the metadata.tsv file can be found in the manual.

8. **Run build_database.py. Detailed instructions on build_database.py can be found in the manual.**

### Some notes about sequence headers:
* Each sequence header (sequence header = Unique ID) within a file must be unique. Sequence
headers must be the same across all files for a taxon and must be the same as the Unique ID
provided in the metadata.tsv file for the taxon. This is true for both ortholog and paralog
fasta files.
* Sequence headers cannot contain underscores “_”, at symbols “@” white spaces, or double
dots “..”. This is true for both ortholog and paralog fasta files.
* If you provided separate sequence files containing paralogs, each sequence header within each
file will have “..p<randomfivedigitnumber>” appended to the end by build_database.py (the
python script that will prepare the custom database for use in PhyloFisher).