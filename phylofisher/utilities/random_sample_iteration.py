"""
Useage:
python makeIQ.py
"""
import os
import sys
import glob
import random
totalgenes = int(sys.argv[1]) #### How many total genes do you have?
iterations = int(sys.argv[2]) #### How many iterations would you like to run
tenpercent = totalgenes*0.1
twentypercent = totalgenes*0.2
thirtypercent = totalgenes*0.3
fortypercent = totalgenes*0.4
fiftypercent = totalgenes*0.5
sixtypercent = totalgenes*0.6
seventypercent = totalgenes*0.7
eightypercent = totalgenes*0.8
ninetypercent = totalgenes*0.9
genes = [fname [:-9] for fname in glob.glob('*.bmge.fas')]
cd=os.getcwd()
for n in range(iterations):
os.system("rm -r twenty/ forty/ sixty/ eighty/")
twenty = random.sample(genes, int(twentypercent))
forty = random.sample(genes, int(fortypercent))
sixty = random.sample(genes, int(sixtypercent))
eighty = random.sample(genes, int(eightypercent))
os.system("mkdir twenty")
for i in twenty:
os.system("cp %s.bmge.fas twenty" %(i))
os.system("cp alvert_septable.py seqtools.py twenty")
os.chdir('twenty')
os.system("python alvert_septable.py -c bmge.fas Bordor.20.%s.dat" %(n))
os.system("cp *.dat* ../")
os.chdir(cd)
os.system("mkdir forty")
for i in forty:
os.system("cp %s.bmge.fas forty" %(i))
os.system("cp alvert_septable.py seqtools.py forty")
os.chdir('forty')
os.system("python alvert_septable.py -c bmge.fas Bordor.40.%s.dat" %(n))
os.system("cp *.dat* ../")
os.chdir(cd)
os.system("mkdir sixty")
for i in sixty:
os.system("cp %s.bmge.fas sixty" %(i))
os.system("cp alvert_septable.py seqtools.py sixty")
os.chdir('sixty')
os.system("python alvert_septable.py -c bmge.fas Bordor.60.%s.dat" %(n))
os.system("cp *.dat* ../")
os.chdir(cd)
os.system("mkdir eighty")
for i in eighty:
os.system("cp %s.bmge.fas eighty" %(i))
os.system("cp alvert_septable.py seqtools.py eighty")
os.chdir('eighty')
os.system("python alvert_septable.py -c bmge.fas Bordor.80.%s.dat" %(n))
os.system("cp *.dat* ../")