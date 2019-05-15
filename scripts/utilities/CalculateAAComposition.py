import os
import sys
import re  
from decimal import *
### AGPST, C, FWY, HRK, MILV, NDEQ === DAYHOFF classes 0-5

infile = sys.argv[1]
out = sys.argv[2]


infile = open(infile, "r")
out = open(out, "w")

lines = infile.readlines()

out.write("taxon\tA\tG\tP\tS\tT\tC\tF\tW\tY\tH\tR\tK\tM\tI\tL\tV\tN\tD\tE\tQ\n")

for line in lines[1:]:
	line =line.strip()
	taxon= line.split()[0]
	sentence=line.split()[1]
	length=len(sentence)
	gaps=sentence.count('-')
	sites = length-gaps
	print sites
	out.write(taxon)
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('A'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('G'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('P'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('S'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('T'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('C'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('F'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('W'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('Y'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('H'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('R'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('K'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('M'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('I'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('L'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('V'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('N'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('D'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('E'))/Decimal(sites))))
	out.write("\t")
	out.write(str(float(Decimal(sentence.count('Q'))/Decimal(sites))))
	out.write("\n")
out.close()


"""
mydata = read.table("/Users/mbrown/Documents/MyDocuments-Aine/Molecular-AINE/Bordor/Bordor.351.64.3-7-2016.Roger/CalculateAAComposition/Bordor.351.64.3-7-2016.dat.AAcomposition.tsv",header=T,row.names = 1)
d <- dist(as.matrix(mydata))
hc <- hclust(d)  
plot(hc)

"""
	
	