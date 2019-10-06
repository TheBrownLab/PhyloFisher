# Example for Divvier
This folder has an example MSA and divvier output to demonstrate how to run the program

### BB20005.subset.fas
An input file that contains a subset of BB20005 from BAliBASE4.0
### BB20005.subset.linsi.fas
An MSA file produced by running linsi from MAFFT on the file above.
```
linsi BB20005.subset.fas > BB20005.subset.linsi.fas
```
### BB20005.subset.linsi.divvier.fas
A divvier file created by running divvier with the options on the homepage
```
./divvier -mincol 4 -divvygap BB20005.subset.linsi.fas
```
### BB20005.subset.linsi.partial.fas
A divvier file created by running divvier with the partial filtering options on the homepage
```
./divvier -partial -mincol 4 -divvygap BB20005.subset.linsi.fas
```

### Note
Running any of the above divvier commands will also create a file
```
BB20005.subset.linsi.fas.PP
```
This file contains the posterior probabilities (PP) computed during the run. This is often the most time consuming step when running divvier, so this file allows you to run it and try different options without having to do that step again. If you do not intend on rerunning and testing divvier then it's safe to remove this file. 
