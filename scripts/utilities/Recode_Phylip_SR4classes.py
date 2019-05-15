import os
import sys
import re  

### AGPST, C, FWY, HRK, MILV, NDEQ === DAYHOFF classes 0-5


"""
A,G,N,P,S,T = A; 
C,H,W,Y = C; 
D,E,K,Q,R = G; 
F,I,L,M,V = T 
"""



infile = sys.argv[1]
out = sys.argv[2]


infile = open(infile, "r")
out = open(out, "w")

lines = infile.readlines()
out.write(lines[0])
for line in lines[1:]:
	#print line
	line = line.strip()
	name = line.split()[0]
	seq = line.split()[1]
	#print seq
	seq = seq.replace("A","A")
	seq = seq.replace("G","A")
	seq = seq.replace("N","A")
	seq = seq.replace("P","A")
	seq = seq.replace("S","A")
	seq = seq.replace("T","A")
	
	seq = seq.replace("C","C")
	seq = seq.replace("H","C")
	seq = seq.replace("W","C")
	seq = seq.replace("Y","C")
	
	seq = seq.replace("D","G")
	seq = seq.replace("E","G")
	seq = seq.replace("K","G")
	seq = seq.replace("Q","G")
	seq = seq.replace("R","G")

	seq = seq.replace("F","T")
	seq = seq.replace("I","T")
	seq = seq.replace("L","T")
	seq = seq.replace("M","T")
	seq = seq.replace("V","T")
	
	seq = seq.replace("X","-")
	out.write(name + " " + seq + "\n")

out.close()
	

