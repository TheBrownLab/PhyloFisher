# PREQUAL - A pre-alignment quality filter
Written by Simon Whelan in collaboration with Fabien Burki and Iker Irisarri

A program to identify and mask regions with non-homologous adjacent characters in FASTA files. The following is just the -h all option with some more description. This problem uses extensive code from [Zorro](https://phylogenomics.me/software/zorro/) and [PaHMM-Tree](https://github.com/marbogusz/paHMM-Tree) and is distributed under a GPL v3.0. Please see the manual for more detailed information about functionality and settings.

### Options affecting the core region and filtering:
	-corerun X       	: X number of high posterior sites at beginning and end before 
					a core region is defined [DEFAULT 3]
	-nocore          	: No core region will be defined [DEFAULT uses core region]
	-removeall       	: Remove all sites rather than those outside the core region 
					[DEFAULT sites in core filtered with X]
	-corefilter X    	: The character outputted in the core for the  filtered alignment [DEFAULT is X]

### Options affecting output formats:
	-outsuffix X     	: Output file will be the original name with X as a suffix [DEFAULT .filtered]
	-dosummary X     	: Output summary statistics to file suffixed with X [if no X then DEFAULT .summary]
	-dodetail X      	: Output detailed statistics to file suffixed with X [if no X then DEFAULT .detail]
	-noPP            	: Stop outputting the posterior probability matrix

### Options affecting posterior probabilities and filtering:
	-pptype X [Y]       	: Specify the algorithm used to calculate posterior probabilities
					X = all : for all against all sequence comparisons
					X = closest : for Y closest relatives [DEFAULT; Y = 10]
					X = longest : for comparing the Y longest sequences [Y = 10]
	-filterprop X       	: Filter the sequences so that in total X proportion (range 0.0 - 1.0) 
					of the sequences are maintained.
	-filterthresh X     	: Filter the sequences to the posterior probabilities threshold X 
					(range 0.0 - 1.0). Any amount of sequence can be removed. 
					Recommended option (DEFAULT X = 0.994).
	-filterjoin X       	: Extend filtering over regions of unfiltered sequence less than X [DEFAULT X = 25]
	-nofilterlist X     	: Specify a file X that contains a list of taxa names that will 
					not be filtered. In X one name per line.
	-nofilterword X     	: Specify a file X that contains a list of words and sequence names that contain 
					those words will not be filtered. In X one word per line.

### Usage:
	./prequal [options] input_file
	-h [all] for [full] options

### Typical usage (should do a good job with most sequences):
	 ./prequal input_file

### Citation
PREQUAL: a program for detecting errors in sets of unaligned homologous sequences. Whelan, Irisarri & Burki. In Press. Bioinformatics. <https://doi.org/10.1093/bioinformatics/bty448>
