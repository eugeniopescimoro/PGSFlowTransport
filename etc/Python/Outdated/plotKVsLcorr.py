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
from mpl_toolkits.mplot3d import axes3d, Axes3D
#import matplotlib.image as img
from pathlib import Path
from Parse import parseLog

# Initialization ##############################################################
Lcorr = []
Mvel = []
sim = 1 # Number of simulations to analyse 
interval = 1 # interval between simulations

# Draw blank canvas ###########################################################
font = {'size': 24}
plt.rc('font', **font)
fig = plt.figure(figsize=(14, 14))

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

# Plot a 3d representation of correlation lengths and the average Vx ##########
ax = fig.add_subplot(2, 2, 1, projection='3d') 
ax.plot3D(x[:10], y[:10], z[:10], 'ro')
ax.plot3D(x[11:20], y[11:20], z[11:20], 'bo')
ax.plot3D(x[21:sim], y[21:sim], z[21:sim], 'go')
ax.set_xlabel('\n' + 'Lx [m]', linespacing = 2, fontdict = font)
ax.set_ylabel('\n' + 'Ly [m]', linespacing = 2, fontdict = font)
ax.set_zlabel('\n' + 'Lz [m]', linespacing = 2, fontdict = font)
ax.set_title('Correlation lengths', fontdict = font)
ax.tick_params(labelsize = 24)
ax = fig.add_subplot(2, 2, 2)
ax.plot(list([i for i in range(sim)]), list(map(itemgetter(0), Mvel[:sim])), 'ro')
#ax.plot(list([i for i in range(11, 20)]), list(map(itemgetter(0), Mvel[11:20])), 'bo')
#ax.plot(list([i for i in range(21, sim)]), list(map(itemgetter(0), Mvel[21:sim])), 'go')
ax.set_xlabel('Simulation number [-]', fontdict=font)
ax.set_ylabel('Average Vx [m/s]', fontdict = font)
ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
ax.tick_params(labelsize = 24)

# Plot the K and velocity field taken from png images #########################
#ax = fig.add_subplot(2, 2, 1)
#Kfield = img.imread('/home/pmxep5/Pictures/veLcorr/TS15/TS15K.png')
#plt.axis('off')
#plt.imshow(Kfield)
#ax.text(0.05, 0.95, "A", transform=ax.transAxes, ha="left", va="top", fontsize=24)
#ax = fig.add_subplot(2, 2, 2)
#Kfield = img.imread('/home/pmxep5/Pictures/veLcorr/TS15/TS15p.png')
#plt.axis('off')
#plt.imshow(Kfield)
#ax.text(0.05, 0.95, "B", transform=ax.transAxes, ha="left", va="top", fontsize=24)

# Plot average Vy and Vz against the simulation number ########################
ax = fig.add_subplot(2, 2, 3)
ax.plot(list([i for i in range(sim)]), list(map(itemgetter(1), Mvel[:sim])), 'ro')
#ax.plot(list([i for i in range(11, 20)]), list(map(itemgetter(1), Mvel[11:20])), 'bo')
#ax.plot(list([i for i in range(21, sim)]), list(map(itemgetter(1), Mvel[21:sim])), 'go')
ax.set_xlabel('Simulation number [-]', fontdict=font)
ax.set_ylabel('Average Vy [m/s]', fontdict = font)
ax.text(0.05, 0.95, "C", transform=ax.transAxes, ha="left", va="top", fontsize=24)
ax.tick_params(labelsize = 24)
ax = fig.add_subplot(2, 2, 4)
ax.plot(list([i for i in range(sim)]), list(map(itemgetter(2), Mvel[:sim])), 'ro')
#ax.plot(list([i for i in range(11, 20)]), list(map(itemgetter(2), Mvel[11:20])), 'bo')
#ax.plot(list([i for i in range(21, sim)]), list(map(itemgetter(2), Mvel[21:sim])), 'go')
ax.set_xlabel('Simulation number [-]', fontdict=font)
ax.set_ylabel('Average Vz [m/s]', fontdict = font)
ax.text(0.05, 0.95, "D", transform=ax.transAxes, ha="left", va="top", fontsize=24)
ax.tick_params(labelsize = 24)

# Adjust visualization, save and show plots ###################################
plt.tight_layout()
fig.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/aveVelCorr.pdf")
#plt.show()
