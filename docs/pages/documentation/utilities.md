---
layout: default
title: Utilities
parent: Documentation
nav_order: 3
permalink: /documentation/utilities
---
# Utilities

* **aa_comp_calculator.py**
    * Calculates amino acid composition and uses euclidean distances to hierarchically
        cluster these data, in order to examine if amino acid composition may bias the groupings that were
        inferred in a phylogenomic tree. See Brown et al., 2018 for an example.
    * Usage: `aa_comp_calculator.py [OPTIONS] -i <input_matrix>`
  
* **astral_runner.py**
    * Generates input files and infers a coalescent-based species tree given a set of single
        ortholog trees and bootstrap trees using ASTRAL-III (Zhang et al., 2018).
    * Usage: `astral_runner.py [OPTIONS] -i <input_directories>`
    
* **backup_restoration.py** 
    * Restores a previous version of the database from the directory backups/.
    * Usage: `backup_restoration.py [OPTIONS] -d <path/to/database/>`
    
* **bipartition_examiner.py** 
    * Calculates the observed occurrences of clades of interest in bootstrap trees.
    * Usage: `bipartition_examiner.py [OPTIONS] -b <input_MLBS_files> -g groups.txt`
    
* **build_database.py** 
    * used to format custom databases for use in the PhyloFisher workflow. This utility
        is also used to rename taxa in either the provided database or a custom database.
      * Usage: `build_database.py [OPTIONS]`
  
* **explore_database.py** 
    * Used to examine the composition of the database using taxonomic terms as search
        queries.
      * Usage: `explore_database.py [OPTIONS]`
  
* **fast_site_remover.py** 
    * The fastest evolving sites are expected to be the most prone to phylogenetic signal
        saturation and systematic model misspecification in phylogenomic analyses. This tool will remove the
        fastest evolving sites within the phylogenomic supermatrix in a stepwise fashion, leading to a user
        defined set of new matrices with these sites removed.
      * Usage: `fast_site_remover.py [OPTIONS] -m <input_matrix>`
  
* **fast_taxa_remover.py** 
    * Removes the fastest evolving taxa, based on branch length. This tool will
        remove the fastest evolving taxa within the phylogenomic supermatrix in a stepwise fashion, leading
        to a user defined set of new matrices with these taxa removed.
      * Usage: `fast_taxa_remover.py [OPTIONS] -m <input_matrix> -t <input_tree>`
  
* **genetic_code_examiner.py**     
    * Checks stop-to-sense and sense-to-sense codon reassignment signal in
        transcriptome/genome data. This script is using the advantage of the phylogenetically broad database
        accompanying our software. The first step is the creation of multiple sequence alignments from all
        fasta files containing manually curated orthologous sequences. This step can take some time but it is
        necessary only for the first run of genetic_code.py. In the next step, tblastn searches with orthologs
        from selected related organism/s in the database are performed against the given transcriptome/genome
        (with 1e-30 default e-value). For all genes, the best scoring hit is investigated for “good quality
        positions”. Such positions are located at least 6 amino acids from the beginning or end of the blast
        alignment and the number of low scoring mismatches (normally denoted by spaces in blast middle line)
        in close proximity to these positions (+- 3 amino acids) is less than 3. Corresponding positions from the
        query are then analyzed in previously prepared multiple sequence alignments and information about
        well-conserved amino acids (in more than 70% of organisms, default) is collected and connected to the
        underlying codon from transcriptome/genome. Codons which show evidence for signals different from
        the standard genetic code signal are then visualized in the form of bar plots (occurrence of conserved
        amino acids). Thanks to the evolutionary well-conserved nature of proteins in our database, realigning
        all sequences again with provided input nucleotide data is not necessary. This script performs well with
        genomic and transcriptomic data. Analysis of one transcriptome/genome should usually take less than
        5 minutes on an average personal computer. It has to be mentioned that alternative genetic code signal
        from multiple sequence alignments is only one way to analyze this phenomenon and tRNAs should be
        investigated as well if possible.
      * Usage: `genetic_code_examiner.py [OPTIONS] -i <input_file> -q <query_file>`
  
* **heterotachy.py** 
    * Within-site rate variation (heterotachy) (Lopez et al., 2002) has been shown to cause
        artificial relationships in molecular phylogenetic reconstruction (Inagaki et al., 2004). This tool will
        remove the most heterotacheous sites within a phylogenomic supermatrix in a stepwise fashion, leading
        to a user defined set of new matrices with these sites removed.
    * Usage: `heterotachy.py -t <input_tree> -m <input_matrix> [OPTIONS]`
  
* **mammal_modeler.py** 
    * Generates a MAMMaL site heterogeneous model from a user input tree and
        supermatrix with estimated frequencies for a user defined number of classes. This program creates a
        set of temporary files from the user provided input that MAMMaL is able to handle. The output is
        a heterogeneous model in Nexus format that is usable in IQtree using options (-m LG+ESmodel+G
        -mdef esmodel.nex). This program outputs by default a 61 class mixture model with 60 site frequency
        classes and the overall amino acid frequencies (+F) of the phylogenomic dataset. Also please note
        that when using this program we hard-code the “not using likelihood weighting” option. This is
        because using likelihood weighting will cause issues in the calculation of likelihood weights in sparse
        phylogenomic matrices. The problem is that there may be pairs of sequences that have no sites in
        common. Specifically, the proportion of times an amino acid occurs for the pair of sequences becomes
        NA because the denominator is 0 (calculation is [p_{aa;sj}] in Eqn (5) of the Susko et al., 2018) (Ed
        Susko personal communication).
      * Usage: `mammal_modeler.py [OPTIONS] -i <input_matrix>]`
  
* **purge.py** 
    * This tool is used for deleting taxa and/or taxonomic groups from the database and metadata
    permanently.
      * Usage: `purge.py [OPTIONS] -i <input_directory>`
  
* **random_resampler.py**
    * This tool randomly resamples the gene set into a set of new matrices that
        are subsamples of the super matrix. It constructs supermatrices from randomly sampled genes with
        user defined options such as the confidence interval sampling all genes in a random fashion and the
        percentage of subsampling a user requires per sampling step. This method was used in Brown et al.,
        2018 as an example.
      * Usage: `random_resampler.py [OPTIONS] -i <input_directory>`
  
* **rtc_binner.py**
    * Calculates the relative tree certainty score (RTC) in RAxML Stamatakis, 2014 of each
        single ortholog tree and bins them based on their RTC scoring into top 25%, 50%, and top 75% sets.
        Supermatrices are constructed from these bins of orthologs.
      * Usage: `rtc_binner.py [OPTIONS] -i <input_matrix>`
  
* **SR4_class_recoder.py**
    * To minimize phylogenetic saturation this tool recodes input supermatrix into
        the four-character state scheme of SR4 (Susko and Roger, 2007), based on amino acid classification.
      * Usage: `SR4_class_recoder.py [OPTIONS] -i <input_matrix>`
  
* **taxon_collapser.py**
    * Allows users to combine multiple operational taxonomic units into one single taxon.
        For example if a user has multiple single cell libraries from a taxon or multiple strains of the same
        species (or genus etc.), a user may decide to collapse all these strains/libraries into a single taxon.
      * Usage: `taxon_collapser.py -i <to_collapse.tsv>`
  
***Use `<script>.py --help`  to see all options and their desctriptions.
 