#import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#import os
import re
for i in range(1, 2):
    # os.chdir('TS%d' % i) # It changes the work directory but it is safer to work with absolute paths
    mvel = []
    conc = []
    time = []
    dimX = []
    dimY = []
    dimZ = []
    with open("/home/pmxep5/OpenFOAM/others/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Herten/punctualInj/system/blockMeshDict") as BMD:       
        for line in BMD:
            if "L" in line:
                dimX.append(float(re.split(' |\;', line)[1]))
                break
        for line in BMD:
            if "D" in line:
                dimY.append(float(re.split(' |\;', line)[1]))
                break
        for line in BMD:
            if "H" in line:
                dimZ.append(float(re.split(' |\;', line)[1]))
                break
    with open("/home/pmxep5/OpenFOAM/others/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Herten/punctualInj/log") as log:
        for _ in range(441): # It skips the first n lines to avoid some misleading "Time =" occurencies
            next(log)
        for line in log: # It greps all the relevant data from the log file
            if "Mean vel =" in line:
                mvel.append([float(i) for i in re.split(' |\(|\)', line)[slice(4, 7)]])
            if "Flux out =" in line:
                conc.append(float(line.split()[3])*float(dimY[0])*float(dimZ[0])*mvel[0][0]) # The flux weighted average given by "Flux Out" it is multiplied by the outlet boundary area and the average longitudinal velocity
            if "Time =" in line:
                time.append(float(line.split()[2]))
    #cBoolean = np.logical_and(np.array(conc)>1e-12, np.array(conc)<1) 
    #c = [val for i, val in enumerate(conc) if cBoolean[i]]
    #t = [val for i, val in enumerate(time) if cBoolean[i]]
    #Cmax = max(conc)
    C0 = 1e-6 # This is the concentration term Su provided in the fvOptions file and its unit is [1/s] 
    Cnorm = [x/C0 for x in conc]
    Tadv = float(dimX[0])/mvel[0][0]
    Tnorm = [val/Tadv for i, val in enumerate(time)]# if cBoolean[i]]

for i in range(0, len(Tnorm), 100):
    T = Tnorm[0:i]
    C = Cnorm[0:i]
    plt.plot(T, C, color='blue', linewidth='6')
    plt.xlim([0, max(Tnorm)])
    plt.ylim([0, 1])
    plt.rcParams.update({'font.size': 18})
    plt.xlabel('$t/t* [-]$')
    plt.ylabel('$c/c_{0} [-]$')
    plt.grid(b=True, which='both', axis='both')
    plt.tight_layout()
    plt.savefig('/home/pmxep5/OneDrive/Nottingham/Results/Images/punctualInj/cTiffs/BTC_time%d.tiff'%(i/100))
