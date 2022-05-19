#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: A BTC for each time step is plotted and saved in tiff format
"""

import os
import matplotlib.pyplot as plt
from pathlib import Path
from bashParse import bashParseLog, parseLog, parseConstants, parseInitialConditions, parseSetFieldsDict #--> This option requires the OpenFOAM log to be written from latest version of adaptiveScalarTransportFoam solver
import numpy as np

sim = 1 # Number of simulations to analyse 
FS = 4 # Number of the First Simulation to analyse
s = 1 # Moving average window size -> smoothness factor
c = []
cl = []
dd = []
t = []
m = []
mvel = []

simPath = ['stopConcAdapTmstp/scat_5-lowContrast/TS4']
homeFolderPath = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/', simPath[0]))
saveFolderPath = Path('/home/pmxep5/OneDrive/Nottingham/Results/Images/uniformInj/BTC/')

# Parse #######################################################################
bashParseLog(sim, FS, homeFolderPath) # Import bashParse.py to use bashParseLog
# parseLog function parses the log file from OpenFOAM and stores the relevant data in different lists
kvol, kval = parseLog(homeFolderPath, s, cl, dd, mvel, c, t, m)
tList = t[0].tolist()

Tadv = float(dd[0][0])/mvel[0][0]
Tnorm = np.array([val/Tadv for i, val in enumerate(t)])# if cBoolean[i]]

processorPath = '/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_5-lowContrast/TS4/processor0'
dirlist = [item for item in os.listdir(processorPath) if os.path.isdir(os.path.join(processorPath, item))]
dirList = []
for x in dirlist:
    try:
        float(x)
        dirList.append(x)
    except ValueError:
        pass
dirList = [float(x) for x in dirList]
dirList.sort()

font = {'size': 50}
plt.rc('font', **font)
plt.figure(figsize=(14, 9))
plt.xlim([0, max(Tnorm[0])])
plt.ylim([0, max(c[0])])
plt.xlabel('$T [-]$')
plt.ylabel('$\overline{c} [-]$')
plt.grid(b=True, which='both', axis='both')
plumesImageIndx = [tList.index(x) for x in tList if x in dirList]
for i in range(len(plumesImageIndx)):
    T = Tnorm[0][0:plumesImageIndx[i]]
    C = c[0][0:plumesImageIndx[i]]
    plt.plot(T, C, color='blue', linewidth='6')
    plt.tight_layout()
    plt.savefig(os.path.join(saveFolderPath, 'BTC_time%d.tiff'%i))