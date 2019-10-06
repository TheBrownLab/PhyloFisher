## Divvier: a program for removing MSA uncertainty by Simon Whelan

Divvier is now submitted for publication. The following options are available in divvier

Clustering options:
```
	-divvy       : do standard divvying (DEFAULT)
	-partial     : do partial filtering by testing removal of individual characters
	-thresh X    : set the threshold for divvying to X (DEFAULT = 0.801)
```
Approximation options: 
```
	-approx X    : minimum number of characters tested in a split during divvying (DEFAULT X = 10)
	-checksplits : go through sequence and ensure there's a pair for every split. Can be slow
	-HMMapprox   : Do the pairHMM bounding approximation (DEFAULT)
	-HMMexact    : Do the full pairHMM and ignore bounding
```
Output options: 
```
	-mincol X    : Minimum number of characters in a column to output when divvying/filtering (DEFAULT X = 2)
	-divvygap    : Output a gap instead of the static * character so divvied MSAs can be used in phylogeny
```

#### Recommended usage

There are two main modes for running divvier on an alignment file: full divvying and partial filtering. Full divvying is the default option and will aim to retain as much information in the MSA as possible. For immediate phylogenetic use we suggest the following that will require at least 4 character per column and produce output to myfile.divvy.fas that will work as input into standard phylogeny programs:
```
./divvier -mincol 4 -divvygap myfile.fas
```
Partial filtering can be done by adding a single option to this command line:
```
./divvier –partial -mincol 4 -divvygap myfile.fas
```


