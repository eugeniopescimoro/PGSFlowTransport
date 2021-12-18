#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Plot average velocities stored in output log file from OpenFOAM simpleDarcyFoam solver against simulation number
"""

# Import section ############################################################## 
from operator import itemgetter
import matplotlib.pyplot as plt
from pathlib import Path
from Parse import parseLog

# Initialization ##############################################################
Lcorr = []
Mvel = []
sim = 1 # Number of simulations to analyse 
interval = 1 # interval between simulations

# Loop through simulation folders and parse files (setRandomFieldDict, log and corrLengths.mat)
for i in range(1, sim+1, interval): # Range of simulations which should be plotted
# Paths 
    luin67272TS = Path('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_2/TS%d' % i)
    #luin67272TS = Path('/home/pmxep5/OpenFOAM/others/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Python/veLcorrTransport/TS%d' % i) 
    gercphdTS = Path('/data/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Python/veLcorrTransport/TS%d' % i)

# Parse section ###############################################################
# The parseLog function parses the log file from OpenFOAM and stores the relevant data in different lists
    cl, dd, mass, mvel, conc, time = parseLog(luin67272TS)
    Lcorr.append(cl[0])
    Mvel.append(mvel)

x = list(map(itemgetter(0), Lcorr))
y = list(map(itemgetter(1), Lcorr))
z = list(map(itemgetter(2), Lcorr))

# Plot average Vx, Vy and Vz against the simulation number ####################
font = {'size': 24}
plt.rc('font', **font)
fig1 = plt.figure(figsize=(14, 14))
ax = fig1.add_subplot(1, 1, 1)
ax.plot(list(map(itemgetter(0), Lcorr[:10])), list(map(itemgetter(0), Mvel[:10])), 'ro')
ax.set_xlabel('Correlation length [m]')
ax.set_ylabel('Average Vx [m/s]')
ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
ax.set_ylim([0.5e-5, 1.1e-5])
plt.tight_layout()
fig1.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/aveVxNonPeriodic.pdf")
#plt.show()
fig2 = plt.figure(figsize=(14, 14))
ax = fig2.add_subplot(1, 1, 1)
ax.plot(list(map(itemgetter(1), Lcorr[10:20])), list(map(itemgetter(1), Mvel[10:20])), 'bo')
ax.set_xlabel('Correlation length [m]')
ax.set_ylabel('Average Vy [m/s]')
ax.set_ylim([-5.5e-7, 3.5e-7])
plt.tight_layout()
fig2.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/aveVyNonPeriodic.pdf")
#plt.show()
fig3 = plt.figure(figsize=(14, 14))
ax = fig3.add_subplot(1, 1, 1)
ax.plot(list(map(itemgetter(2), Lcorr[20:sim])), list(map(itemgetter(2), Mvel[20:sim])), 'go')
ax.set_xlabel('Correlation length [m]')
ax.set_ylabel('Average Vz [m/s]')
ax.set_ylim([-5.5e-7, 3.5e-7])
plt.tight_layout()
fig3.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/aveVzNonPeriodic.pdf")
#plt.show()