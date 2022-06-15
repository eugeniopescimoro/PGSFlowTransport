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
from Parse import parseLog, parseConstants, parseInitialConditions, parseSetFieldsDict
from pathlib import Path
from Process import processConc

# Initialization ##############################################################
sim = 1 # Number of simulations to analyse 
interval = 1 # Interval between increasing simulations
s = 100 # Smoothness plotting factor: it reads one value every "s"

# Draw blank canvas, grids and legend #########################################
font = {'size': 24}
plt.rc('font', **font)
fig1, (f1ax1, f1ax2, f1ax3)  = plt.subplots(1, 3, figsize=(14, 9))
fig2, (f2ax1, f2ax2, f2ax3) = plt.subplots(1, 3, figsize=(14, 9))
fig3, (f3ax1, f3ax2, f3ax3) = plt.subplots(1, 3, figsize=(14, 9))
ax = plt.gca() # Get Current Axis to access the line colours with ax._get_lines.prop_cycler
for axs in f1ax1, f1ax2, f1ax3, f2ax1, f2ax2, f2ax3, f3ax1, f3ax2, f3ax3:
    axs.grid(True)
f1ax1.set_xlabel('$t/t* [-]$')
f1ax1.set_ylabel('$c [-]$')
f2ax1.set_xlabel('$t/t* [-]$')
f2ax1.set_ylabel('$dc/dt [-]$')
f3ax1.set_xlabel('$t/t* [-]$')
f3ax1.set_ylabel('$1-c [-]$')

# Loop over the simulations ################################################### 
for i in range(1, sim+1, interval):
# Paths 
    homeFolderPath = Path('/data/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_3/TS%d' % i)
    # homeFolderPath = Path('/data/PGSFlowTransport/tutorials/RESULTS/Herten/Herten9')
    # homeFolderPath = Path('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_3/TS%d' % i)
# Parse section ###############################################################
# parseLog function parses the log file from OpenFOAM and stores the relevant data in different lists
    cl, dd, mass, mvel, conc, time = parseLog(homeFolderPath, s)
# parseTransportProperties function parses the transportProperties OpenFOAM input file and stores the hydraulic dispersion coefficient value
    D, mu, rho, g = parseConstants(homeFolderPath)
# parseFieldsDict function parses the setFieldsDict dictionary and outputs the X distance at which the 
# center of the injection volume is placed (not needed if the injection is on the inlet boundary) 
    topoPath = os.path.join(homeFolderPath, "system/topoSetDict")
    if os.path.isfile(topoPath):
        Xbox = parseSetFieldsDict(homeFolderPath)
    else:
        Xbox = 0
# parseInitialConditions function parses the initial pressure conditions
    deltaPx = parseInitialConditions(homeFolderPath)
    
# Processing section ##########################################################
    y, mu1, mu2, lam, c, t, dC, dY = processConc(homeFolderPath, dd, mvel, conc, time, Xbox, s)
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
    if not cl:
        cl = [[0, 0, 0]] # In case the permeability field is generated out of OpenFOAM and it has unknown correlation length     
    medPeX = cl[0][0]*mvel[0]/D[0]
    medPeY = cl[0][1]*mvel[1]/D[0]
    medPeZ = cl[0][2]*mvel[2]/D[0]
    # Concentration vector norm
    normC = np.linalg.norm(c)
    normY = np.linalg.norm(y)
    # Effective permeatbility
    Kxeff = -mvel[0]*mu[0]/(deltaPx-rho[0]*g[2])
    
# Print section ###############################################################
    # From log    
    with open("log") as stats: # It re-opens the log file to print the field statistics 
        for line in stats:
            if "Gaussian Random Field generated" in line:
                print("============= SIMULATION %d =============\n\nGEOSTATISTICAL METRICS" % i)
                while set(line.split()).isdisjoint(set(["End"])):
                    if "Statistics log" in line:
                        print(line, end = ''),
                        line = stats.readline()
                    print(line, end = '')
                    line = stats.readline()   
    print("Effective longitudinal permeability Kx,eff = %.15f \n" % Kxeff)
    print("FLOW METRICS")
    print("Mean Vx Vy Vz = %.9f %.9f %.9f" % (mvel[0], mvel[1], mvel[2]))
    print("Estimated mean Vx = %.9f \n" % (v))  
    print("Péclet: \n  Macro = (%.2f %.2f %.2f) \n  Meso = (%.2f %.2f %.2f) \n" % (macroPeX, macroPeY, macroPeZ, medPeX, medPeY, medPeZ))
    print("TRANSPORT METRICS")
    print("mu_1 = %f \nmu_2 = %f \nInv Gau lambda = %f \nalpha_estimatedMeanVx = %f \nalpha_meanVx = %f \n" % (mu1, mu2, lam, alphaE, alpha))                 
    print("Dispersion: \n  Molecular = %.9f \n  Estimated mechanical = %f \n  Estimated hydrodynamic = %f \n" % (D[0], MD, HD))
    print("Computed ||c|| = %f \nEstimated ||c|| = %f \n" % (normC, normY))
    
# Plot section ################################################################
    color=next(ax._get_lines.prop_cycler)['color']   
    
    ax.set_ylim([min(c),max(c)])
    
    f1ax1.plot(t, c, color=color, label='TS%d experimental' % i)
    f1ax2.semilogy(t, c, color=color)
    f1ax3.loglog(t, c, color=color)
    f1ax1.plot(t, y, '--', color=color, label='TS%d cum inv Gau' % i)
    f1ax2.semilogy(t, [i for i in y], '--', color=color)
    f1ax3.loglog(t, [i for i in y], '--', color=color)

    f2ax1.plot(t[:-s], dC, color=color)
    f2ax2.semilogy(t[:-s], dC, color=color)
    f2ax3.loglog(t[:-s], dC, color=color)  
    f2ax1.plot(t[:-s], dY, '--', color=color)
    f2ax2.semilogy(t[:-s], dY, '--', color=color)
    f2ax3.loglog(t[:-s], dY, '--', color=color)  
    
    ax.set_ylim([1-max(c),1-min(c)])
    
    f3ax1.plot(t, [1-i for i in c], color=color)
    f3ax2.semilogy(t, [1-i for i in c], color=color)
    f3ax3.loglog(t, [1-i for i in c], color=color)
    f3ax1.plot(t, 1-y, '--', color=color)
    f3ax2.semilogy(t, [1-i for i in y], '--', color=color)
    f3ax3.loglog(t, [1-i for i in y], '--', color=color)
    
    plt.tight_layout()
    
# Plot c, (1-c) and (dc/dt) on 3 figures, each with linear, semilog and log log

handles, labels = f1ax1.get_legend_handles_labels()
fig1.legend(handles, labels, prop=font)
fig2.legend(handles, labels, prop=font)
fig3.legend(handles, labels, prop=font)
# Adjust visualization, save and show plots
os.makedirs(os.path.join(homeFolderPath, "../images"), exist_ok = True)
fig1.savefig(os.path.join(homeFolderPath, "../images/BTCmomentC.pdf"))
fig2.savefig(os.path.join(homeFolderPath, "../images/BTCmomentDCDT.pdf"))
fig3.savefig(os.path.join(homeFolderPath, "../images/BTCmoment1-C.pdf"))
#fig1.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/BTCmoment.pdf")
plt.show()

###############################################################################
# The scipy.stats.invgauss.cdf is a standard scipy method to compute an Inverse Gaussian function
# however it only accepts one parameter (mu) because it assumes lambda = 1 -> we need 2 parameters 
#import scipy
#anaCdf = scipy.stats.invgauss.cdf(t, mu1)
#axs1[0, 0].plot(t, anaCdf)
#aa = scipy.stats.invgauss.fit(dC)
###############################################################################
