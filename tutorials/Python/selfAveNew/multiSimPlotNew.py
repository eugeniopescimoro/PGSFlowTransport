import numpy as np
from numpy import diff
import matplotlib.pyplot as plt
import os
import re
fig, axs = plt.subplots(2, 2)
for i in range(1, 11):
    os.chdir('/home/pmxep5/OpenFOAM/others/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Python/selfAveNew/TS%d' % i)
    dd = []
    cl = []
    mass = []
    mvel = []
    conc = []
    time = []
    lengthX = []
    with open("log") as log:
        for _ in range(225): # It skips the first 225 to avoid some misleading "Time =" occurencies
            next(log)
        for line in log: # It greps all the relevant data from the log file
            if "Total mass =" in line:
                mass.append(float(line.split()[3]))
            if "Mean vel =" in line:
                mvel.append([float(i) for i in re.split(' |\(|\)', line)[slice(4, 7)]])
            if "Flux out =" in line:
                conc.append(float(line.split()[3]))
            if "Time =" in line:
                time.append(float(line.split()[2]))
    cBoolean = np.logical_and(np.array(conc)>1e-12, np.array(conc)<1) 
    c = [val for i, val in enumerate(conc) if cBoolean[i]]
#    with open("system/blockMeshDict") as BMD:       
#        for line in BMD:
#            if "L" in line:
#                lengthX.append(float(re.split(' |\;', line)[1]))
#                break
#    Tadv = float(lengthX[0])/mvel[0][0]
    Tadv = 2/mvel[0][0]
    t = [val/Tadv for i, val in enumerate(time) if cBoolean[i]]
    dC = diff(c)/diff(t)
    axs[0, 0].plot(t, c)
    axs[0, 0].set_xlabel('$t/t* [-]$')
    axs[0, 0].set_ylabel('$c [-]$')
    axs[0, 1].plot(t[1:], [val for i, val in enumerate(dC)])
    axs[0, 1].set_xlabel('$t/t* [-]$')
    axs[0, 1].set_ylabel('$dc/dt [-]$')
    axs[1, 0].loglog(t, [val for i, val in enumerate(c)])
    axs[1, 0].set_xlabel('$t/t* [-]$')
    axs[1, 0].set_ylabel('$c [-]$')
    axs[1, 1].semilogy(t, [1-val for i, val in enumerate(c)], label='TS%d' % i)    
    axs[1, 1].set_xlabel('$t/t* [-]$')
    axs[1, 1].set_ylabel('$1-c [-]$')
    handles, labels = axs[1, 1].get_legend_handles_labels()
    with open("log") as stats: # It re-opens the log file to print the field statistics 
        for line in stats:
            if "boundingBox:" in line:
                dd.append([float(i) for i in re.split(' |\(|\)', line)[slice(9, 12)]])
            if "Start computing" in line:
                print("=============SIMULATION %d =============\n\n" % i)
                while set(line.split()).isdisjoint(set(["End"])):
                    if "Statistics log" in line:
                        print(line, end = ''),
                        line = stats.readline()
                        cl.append([float(i) for i in re.split(' |\(|\)', line)[slice(2, 5)]])            
                        macroPeX = dd[0][0]*mvel[0][0]/1e-10
                        macroPeY = dd[0][1]*mvel[0][1]/1e-10
                        macroPeZ = dd[0][2]*mvel[0][2]/1e-10
                        medPeX = cl[0][0]*mvel[0][0]/1e-10
                        medPeY = cl[0][1]*mvel[0][1]/1e-10
                        medPeZ = cl[0][2]*mvel[0][2]/1e-10
                        print("P??clet: \n  Macro = (%.2f %.2f %.2f) \n  Meso = (%.2f %.2f %.2f) \n" % (macroPeX, macroPeY, macroPeZ, medPeX, medPeY, medPeZ))                
                    print(line, end = '')
                    line = stats.readline()
    os.chdir('..')
plt.rc('font', size=12)
fig.legend(handles, labels, loc = 'lower center', ncol=10, prop={'size': 6})
axs[0, 0].grid()
axs[0, 1].grid()
axs[1, 0].grid(True, which="both")
axs[1, 1].grid()
plt.tight_layout()
fig.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/selfAveNew/multiSimPlotNew.pdf")
