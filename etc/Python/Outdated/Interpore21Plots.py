#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Plot experimental BTCs from OpenFOAM log and estimate the correspondent cumulative inverse Gaussian function
"""

# Import section ############################################################## 
import numpy as np
from numpy import diff
import matplotlib.pyplot as plt
import os
from Parse import parseLog, parseTransportProperties, parseSetFieldsDict
from pathlib import Path
from Process import processConc

# Draw blank canvas, grids and legend #########################################
font = {'size': 12}
plt.rc('font', **font)
fig1, axs1 = plt.subplots(2, 1, figsize=(14, 9))
ax1 = plt.gca() # Get Current Axis to access the line colours with ax1._get_lines.prop_cycler
axs1[0].grid()
axs1[1].grid()
# axs1[1, 0].grid()
# axs1[1, 1].grid()
# axs1[0, 0].set_xlabel('$t/t* [-]$')
# axs1[0, 0].set_ylabel('$c [-]$')
axs1[0].set_xlabel('$t/t* [-]$')
axs1[0].set_ylabel('$dc/dt [-]$')
# axs1[1, 0].set_xlabel('$t/t* [-]$')
# axs1[1, 0].set_ylabel('$c [-]$')
axs1[1].set_xlabel('$t/t* [-]$')
axs1[1].set_ylabel('$1-c [-]$')

# Loop over the simulations ################################################### 
for i in range(1, 29, 4):
# Paths 
    luin67272TS = Path('/home/pmxep5/OpenFOAM/others/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Python/veLcorrTransport/TS%d' % i) 
    gercphdTS = Path('/data/GMRTFoam/tutorials/simpleDarcyFoam/RESULTS/Python/veLcorrTransport/TS%d' % i)
    luin67272Py = Path('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/etc/Python')
    gercphdPy = Path('/data/PGSFlowTransport/etc/Python')
# Parse section ###############################################################
# parseLog function parses the log file from OpenFOAM and stores the relevant data in different lists
    cl, dd, mass, mvel, conc, time = parseLog(luin67272TS)
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
    HD = D+MD
    alphaE = MD/v
    alpha = MD/mvel[0]
    # Péclet
    macroPeX = dd[0][0]*mvel[0]/D[0]
    macroPeY = dd[0][1]*mvel[1]/D[0]
    macroPeZ = dd[0][2]*mvel[2]/D[0]
    medPeX = cl[0][0]*mvel[0]/D[0]
    medPeY = cl[0][1]*mvel[1]/D[0]
    medPeZ = cl[0][2]*mvel[2]/D[0]
    # Concentration vector norm
    normC = np.linalg.norm(c)
    normY = np.linalg.norm(y)

# Print section ###############################################################
    # From log    
    with open("log") as stats: # It re-opens the log file to print the field statistics 
        for line in stats:
            if "Start computing" in line:
                print("============= SIMULATION %d =============\n\nGEOSTATISTICAL METRICS" % i)
                while set(line.split()).isdisjoint(set(["End"])):
                    if "Statistics log" in line:
                        print(line, end = ''),
                        line = stats.readline()
                    print(line, end = '')
                    line = stats.readline()                
    print("FLOW METRICS")
    print("Mean Vx Vy Vz = %.9f %.9f %.9f" % (mvel[0], mvel[1], mvel[2]))
    print("Estimated mean Vx = %.9f \n" % (v))  
    print("Péclet: \n  Macro = (%.2f %.2f %.2f) \n  Meso = (%.2f %.2f %.2f) \n" % (macroPeX, macroPeY, macroPeZ, medPeX, medPeY, medPeZ))
    print("TRANSPORT METRICS")
    print("mu_1 = %f \nmu_2 = %f \nInv Gau lambda = %f \nalpha_estimatedMeanVx = %f \nalpha_meanVx = %f \n" % (mu1, mu2, lam, alphaE, alpha))                 
    print("Dispersion: \n  Molecular = %.9f \n  Estimated mechanical = %f \n  Estimated hydrodynamic = %f \n" % (D[0], MD, HD))
    print("Computed ||c|| = %f \nEstimated ||c|| = %f \n" % (normC, normY))
    
# Plot section ################################################################
    color=next(ax1._get_lines.prop_cycler)['color']    
    # axs1[0, 0].plot(t, c, color=color)
    # axs1[0, 1].plot(t[1:], [i for i in dC], color=color)
    # axs1[1, 0].loglog(t, [i for i in c], color=color)
    # axs1[1, 1].semilogy(t, [1-i for i in c], label='TS%d' % i, color=color)
    # axs1[0, 0].plot(t, y, '--', color = color)
    # axs1[1, 1].semilogy(t, [1-i for i in y], '--', color=color)
    # axs1[1, 0].loglog(t, [i for i in y], '--', color = color)

    axs1[0].plot(t[1:], [i for i in dC], color=color)
    axs1[1].semilogy(t, [1-i for i in c], label='TS%d' % i, color=color)

handles, labels = axs1[1].get_legend_handles_labels()
fig1.legend(handles, labels, loc = 'lower center', ncol=10, prop={'size': 6})
# Adjust visualization, save and show plots
plt.tight_layout()
#fig1.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/BTCmoment.pdf")
fig1.savefig("/home/pmxep5/OneDrive/shared/BTCmoment.jpg")
#plt.show()

###############################################################################
# The scipy.stats.invgauss.cdf is a standard scipy method to compute an Inverse Gaussian function
# however it only accepts one parameter (mu) because it assumes lambda = 1 -> we need 2 parameters 
#import scipy
#anaCdf = scipy.stats.invgauss.cdf(t, mu1)
#axs1[0, 0].plot(t, anaCdf)
#aa = scipy.stats.invgauss.fit(dC)
###############################################################################
