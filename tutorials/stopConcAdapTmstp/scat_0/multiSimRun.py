#!/usr/bin/env python3
#import numpy as np
import os
import shutil
import re
import subprocess
import pathlib
n = 28
Lcorr = [
 ["0.1", "0.05", "0.05"],
 ["0.2", "0.05", "0.05"],
 ["0.3", "0.05", "0.05"],
 ["0.4", "0.05", "0.05"],
 ["0.5", "0.05", "0.05"],
 ["0.6", "0.05", "0.05"],
 ["0.7", "0.05", "0.05"],
 ["0.8", "0.05", "0.05"],
 ["0.9", "0.05", "0.05"],
 ["1", "0.05", "0.05"],
 ["1", "0.10", "0.05"],
 ["1", "0.15", "0.05"],
 ["1", "0.20", "0.05"],
 ["1", "0.25", "0.05"],
 ["1", "0.30", "0.05"],
 ["1", "0.35", "0.05"],
 ["1", "0.40", "0.05"],
 ["1", "0.45", "0.05"],
 ["1", "0.50", "0.05"],
 ["1", "0.50", "0.10"],
 ["1", "0.50", "0.15"],
 ["1", "0.50", "0.20"],
 ["1", "0.50", "0.25"],
 ["1", "0.50", "0.30"],
 ["1", "0.50", "0.35"],
 ["1", "0.50", "0.40"],
 ["1", "0.50", "0.45"],
 ["1", "0.50", "0.50"]
 ]
kzAniso = ["1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"]
# Remove old directories
for i in range(1, n+1):
    if os.path.exists('TS%d' % i):
        shutil.rmtree('TS%d' % i)
# Copy 'baseCase' directory n times  
for i in range(1, n+1):
    shutil.copytree("baseCase/", "TS%d" % i)
    os.system("cp TS%d/system/setRandomFieldDict.orig TS%d/system/setRandomFieldDict" % (i, i)) 
    os.system("cp TS%d/0/K.orig TS%d/0/K" % (i, i))
    with open('TS%d/system/setRandomFieldDict' % i, 'r') as srfd:
        filedata = srfd.read()
    filedata = filedata.replace('LcorrX', Lcorr[i-1][0])
    filedata = filedata.replace('LcorrY', Lcorr[i-1][1])
    filedata = filedata.replace('LcorrZ', Lcorr[i-1][2])    
    with open('TS%d/system/setRandomFieldDict' % i, 'w') as srfd:
        srfd.write(filedata)            
    with open("TS%d/0/K" % i, 'r') as pr:
        filedata = pr.read()
    filedata = filedata.replace('Kz', kzAniso[i-1])
    with open("TS%d/0/K" % i, 'w') as pr:
        pr.write(filedata) 
    os.chdir("TS%d" % i)
    current_dir = os.getcwd()
    subprocess.run(os.path.join(current_dir,"Allrun"))
    os.chdir("..")
    # It is better to leave it manual for the ease of parallelization 

#########################################################################
################ ORIGINAL Allproperties EXECUTABLE FILE #################
#os.system("cp system/setRandomFieldDict.orig system/setRandomFieldDict") 
#os.system("cp 0/K.orig 0/K")
#os.system("sed -i "s/LcorrX/$1/g" system/setRandomFieldDict") 
#os.system("sed -i "s/LcorrY/$2/g" system/setRandomFieldDict") 
#os.system("sed -i "s/LcorrY/$3/g" system/setRandomFieldDict") 
#os.system("s/Kz/$4/g" 0/K") 
