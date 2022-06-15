#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Plot mechanical dispersion index (alpha) vs non-dimensional Correlation 
Lengths (non-ergodic case) or vs Péclet (ergodic case) from OpenFOAM log file
"""

# Import section ############################################################## 
import numpy as np
from numpy import diff
import os
from operator import itemgetter
import matplotlib.pyplot as plt
from pathlib import Path
from Parse import parseLog, parseTransportProperties, parseSetFieldsDict
from Process import processConc

# Initialization ##############################################################
MechD = []
EalphaX = [] # Estimated alpha X
alphaX = [] 
Evx = [] # Estimated Vx
Mvx = []
sim = 1 # Number of simulations to analyse 
interval = 1 # interval between simulations
Lcorr = [[] for i in range(sim)] # Needed to perform the dot division
Peclet = [[] for i in range(sim)] # Needed to perform the dot division

# Loop through simulation folders and parse files (setRandomFieldDict, log and corrLengths.mat)
for i in range(1, sim+1): # Range of simulations which should be plotted
# Paths 
    luin67272TS = Path('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_2/TS%d' % i)
    #luin67272TS = Path('/home/pmxep5/OpenFOAM/others/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Python/veLcorrTransport/TS%d' % i) 
    gercphdTS = Path('/data/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Python/veLcorrTransport/TS%d' % i)
    luin67272Py = Path('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/etc/Python')
    gercphdPy = Path('/data/PGSFlowTransport/etc/Pytho')

# Parse section ###############################################################
# The parseLog function parses the log file from OpenFOAM and stores the relevant data in different lists
    cl, dd, mass, mvel, conc, time = parseLog(luin67272TS)    
    for number1, number2 in zip(cl[0], dd[0]): Lcorr[i-1].append(number1 / number2) # dot products Matlab equivalent
    Mvx.append(mvel[0])
# parseTransportProperties function parses the transportProperties OpenFOAM input file and stores the hydraulic dispersion coefficient value
    D = parseTransportProperties(luin67272TS)

# parseFieldsDict function parses the setFieldsDict dictionary and outputs the X distance at which the 
# center of the injection volume is placed (not needed if the injection is on the inlet boundary) 
    topoPath = os.path.join(luin67272TS, "system/topoSetDict")
    if os.path.isfile(topoPath):
        Xbox = parseSetFieldsDict(luin67272TS)
    else:
        Xbox = 0
    
# Processing section ##########################################################
    # Meso-Péclet computed using correlation length as characteristic length
    for n1, n2, n3 in zip(cl[0], mvel, D*len(mvel)): Peclet[i-1].append(n1 * n2 /n3)
    # Process concentration data using processConc function 
    y, mu1, mu2, lam, c, t, dC = processConc(luin67272TS, luin67272Py, dd, mvel, conc, time, Xbox)
    # Mean advective velocity and Macro Dispersion estimation from statistical moments (Yu 1999)
    Evx.append(dd[0][0]/mu1)
    MD = dd[0][0]**2/(2*lam)
    MechD.append(MD)
    EalphaX.append(MechD[i-1]/Evx[i-1])
    alphaX.append(MechD[i-1]/Mvx[i-1])

x = list(map(itemgetter(0), Lcorr))
y = list(map(itemgetter(1), Lcorr))
z = list(map(itemgetter(2), Lcorr))

# Plot alpha against non-dimensional correlation length and Péclet ####################
font = {'size': 24}
plt.rc('font', **font)
fig1 = plt.figure(figsize=(14, 9))
ax = fig1.add_subplot(1, 1, 1)
ax.plot(list(map(itemgetter(0), Lcorr[:sim])), EalphaX[:sim], 'ro')
ax.plot(list(map(itemgetter(0), Lcorr[:sim])), alphaX[:sim], 'bx')
ax.set_xlabel('$Lcorr_x/L_x$ [-]')
ax.set_ylabel('Alpha [m]')
ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
ax.set_ylim([min(min(EalphaX[:sim], alphaX[:sim])), max(max(EalphaX[:sim], alphaX[:sim]))])
plt.tight_layout()
fig1.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/ALPHAvsCL_0-sim.pdf")
#plt.show()
fig2 = plt.figure(figsize=(14, 9))
ax = fig2.add_subplot(1, 1, 1)
ax.plot(list(map(itemgetter(0), Peclet[:sim])), EalphaX[:sim], 'ro')
ax.plot(list(map(itemgetter(0), Peclet[:sim])), alphaX[:sim], 'bx')
ax.set_xlabel('Péclet [-]')
ax.set_ylabel('Alpha [m]')
ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
ax.set_ylim([min(min(EalphaX[:sim], alphaX[:sim])), max(max(EalphaX[:sim], alphaX[:sim]))])
plt.tight_layout()
fig2.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/ALPHAvsPECLET_0-sim.pdf")
#plt.show()