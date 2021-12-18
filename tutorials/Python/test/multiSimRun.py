#import numpy as np
import os
import shutil
import re
n = 11
Lcorr = ["0.2", "0.2", "0.1"]
kzAniso = ["1"]
# Remove old directories
for i in range(1, n):
    if os.path.exists('TS%d' % i):
        shutil.rmtree('TS%d' % i)
# Copy 'testPython' directory n times  
for i in range(1, n):
    shutil.copytree("testPython/", "TS%d" % i)
    os.system("cp TS%d/system/setRandomFieldDict.orig TS%d/system/setRandomFieldDict" % (i, i)) 
    os.system("cp TS%d/0/K.orig TS%d/0/K" % (i, i))
    with open('TS%d/system/setRandomFieldDict' % i, 'r') as srfd:
        filedata = srfd.read()
    filedata = filedata.replace('LcorrX', Lcorr[0])
    filedata = filedata.replace('LcorrY', Lcorr[1])
    filedata = filedata.replace('LcorrZ', Lcorr[2])    
    with open('TS%d/system/setRandomFieldDict' % i, 'w') as srfd:
        srfd.write(filedata)            
    with open("TS%d/0/K" % i, 'r') as pr:
        filedata = pr.read()
    filedata = filedata.replace('Kz', kzAniso[0])
    with open("TS%d/0/K" % i, 'w') as pr:
        pr.write(filedata)
    #os.system('gnome-terminal -- ./Allrun') # it will not work unless input from the simulation folder, e.g. TS*/
    # It is better to leave it manual for the ease of parallelization 

#########################################################################
################ ORIGINAL Allproperties EXECUTABLE FILE #################
#os.system("cp system/setRandomFieldDict.orig system/setRandomFieldDict") 
#os.system("cp 0/K.orig 0/K")
#os.system("sed -i "s/LcorrX/$1/g" system/setRandomFieldDict") 
#os.system("sed -i "s/LcorrY/$2/g" system/setRandomFieldDict") 
#os.system("sed -i "s/LcorrY/$3/g" system/setRandomFieldDict") 
#os.system("s/Kz/$4/g" 0/K") 



