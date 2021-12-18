#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Plot Macro-Dispersion vs Correlation Lengths from OpenFOAM log file
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
Lcorr = []
Mvel = []
MechD = []
sim = 1 # Number of simulations to analyse 
interval = 1 # interval between simulations

# Loop through simulation folders and parse files (setRandomFieldDict, log and corrLengths.mat)
for i in range(1, sim+1, interval): # Range of simulations which should be plotted
# Paths 
    luin67272TS = Path('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_2/TS%d' % i)
    #luin67272TS = Path('/home/pmxep5/OpenFOAM/others/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Python/veLcorrTransport/TS%d' % i) 
    gercphdTS = Path('/data/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Python/veLcorrTransport/TS%d' % i)
    luin67272Py = Path('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/etc/Python')
    gercphdPy = Path('/data/PGSFlowTransport/etc/Pytho')

# Parse section ###############################################################
# The parseLog function parses the log file from OpenFOAM and stores the relevant data in different lists
    cl, dd, mass, mvel, conc, time = parseLog(luin67272TS)
    Lcorr.append(cl[0])
# parseTransportProperties function parses the transportProperties OpenFOAM input file and stores the hydraulic dispersion coefficient value
    D = parseTransportProperties(luin67272TS)
# parseFieldsDict function parses the setFieldsDict dictionary and outputs the X distance at which the 
# center of the injection volume is placed (not needed if the injection is on the whole face) 
    Xbox = parseSetFieldsDict(luin67272TS)

# Processing section ##########################################################
    y, mu1, mu2, lam, c, t, dC = processConc(luin67272TS, luin67272Py, dd, mvel, conc, time, Xbox)
    # Mean advective velocity and Macro Dispersion estimation from statistical moments (Yu 1999)
    v = dd[0][0]/mu1
    MD = dd[0][0]**2/(2*lam)
    MechD.append(MD)

x = list(map(itemgetter(0), Lcorr))
y = list(map(itemgetter(1), Lcorr))
z = list(map(itemgetter(2), Lcorr))

# Plot average Vx, Vy and Vz against the simulation number ####################
font = {'size': 24}
plt.rc('font', **font)
fig1 = plt.figure(figsize=(14, 9))
ax = fig1.add_subplot(1, 1, 1)
ax.plot(list(map(itemgetter(0), Lcorr[:10])), MechD[:10], 'ro')
ax.set_xlabel('Correlation length X [m]')
ax.set_ylabel('Mechanical dispersion [m^2/s]')
ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
ax.set_ylim([min(MechD[:10]), max(MechD[:10])])
plt.tight_layout()
fig1.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/MDvsCL_0-10.pdf")
#plt.show()
fig2 = plt.figure(figsize=(14, 9))
ax = fig2.add_subplot(1, 1, 1)
ax.plot(list(map(itemgetter(1), Lcorr[10:20])), MechD[10:20], 'bo')
ax.set_xlabel('Correlation length Y [m]')
ax.set_ylabel('Mechanical dispersion [m^2/s]')
ax.set_ylim([min(MechD[10:20]), max(MechD[10:20])])
plt.tight_layout()
fig2.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/MDvsCL_10-20.pdf")
#plt.show()
fig3 = plt.figure(figsize=(14, 9))
ax = fig3.add_subplot(1, 1, 1)
ax.plot(list(map(itemgetter(2), Lcorr[20:sim])), MechD[20:sim], 'go')
ax.set_xlabel('Correlation length Z [m]')
ax.set_ylabel('Mechanical dispersion [m^2/s]')
ax.set_ylim(min(MechD[20:sim]), max(MechD[20:sim]))
plt.tight_layout()
fig3.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/MDvsCL_20-30.pdf")
#plt.show()

# fig4 = plt.figure(figsize=(14, 9))
# ax2 = fig2.add_subplot(1, 1, 1)
# ax2.plot(cl[0][0], MD, 'o')