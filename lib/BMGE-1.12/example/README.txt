### Entropy-based character trimming ###

Each of the two files prmA.faa and prmA.fna contains a multiple sequence alignment in FASTA format. 
The file  prmA.faa was built from  cyanobacterial prmA amino acid sequences,  and the file prmA.fna 
corresponds to their nucleotide sequences.  The other files prmA.* were obtained with the following 
command lines:

  java -jar BMGE.jar -i prmA.faa -t AA -m BLOSUM95 -o prmA.blosum95.phy -oh prmA.blosum95.html

  java -jar BMGE.jar -i prmA.fna -t DNA -m DNAPAM150:2 -o prmA.pam150.phy -oh prmA.pam150.html



### Stationary-based character trimming ###

The file Rokas.2003.phy is  the supermatrix of character built by Rokas, Williams, King and Carroll 
(2003; Nature, 425:798-804).  The files Rokas.2003.hm.phy and Rokas.2003.hm.html were obtained with
the following command line:

  java -jar BMGE.jar -i Rokas.2003.phy -t DNA -h 1 -g 1 -w 1 -s FAST -o Rokas.2003.hm.phy -oh Rokas.2003.hm.html





