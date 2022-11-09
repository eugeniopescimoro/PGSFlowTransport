#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Plot experimental BTCs from OpenFOAM log and estimate the correspondent cumulative inverse Gaussian function
"""

# IMPORT SECTION ############################################################## 
import numpy as np
import matplotlib.pyplot as plt
import os
# from bashParse import bashParseLog, parseLog, parseConstants, parseInitialConditions, parseSetFieldsDict #--> This option requires the OpenFOAM log to be written from latest version of adaptiveScalarTransportFoam solver 
from Parse import parseLog, parseConstants, parseInitialConditions, parseSetFieldsDict #--> If selected, bashParse and bashParseLog function needs to be commented out
from pathlib import Path
from Process import processConc
from InvGau import invGaussianCDF, invGauVG
from statistics import mean

# INITIALIZATION ##############################################################
sim = 3 # Number of simulations to analyse 
FS = 1 # Number of the First Simulation to analyse
interval = 1 # Interval between increasing simulations
s = 1 # Moving average window size -> smoothness factor
t = [] # List which sotres dimensional times
tt = [] # List which sotres non-dimensional times
c = []
dc = []
dCnorm = []
m = []
dYYnorm = [[]]*sim
y = [[]]*sim
Y = []
yVG = []
kave = []
cl = []
dd = []
kxeff = []
kr = []
HD = []
ndHD = []
mvel = []
Tadv = []
tLS = []
dcLS = []
dcLSnorm = []
kvol = []
kval = []
n = 1 # Derivative smoothing factor

# LOOP THROUGH THE SIMULATIONS ################################################    
for i in range(0, sim, interval):
# Paths
    # simPath = ['/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/testMeshResolution/2D/uniform05cm', '/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/testMeshResolution/2D/uniform1cm', '/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/testMeshResolution/2D/uniform2cm']
    # simPath = ['stopConcAdapTmstp/scat_6-sameDomain/lowCont_seed100', '/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/testMeshResolution/3D/lowCont_seed100_2cm']
    # simPath = ['/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/testMeshResolution/2D/mesh1mm', '/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/testMeshResolution/2D/mesh05cm', '/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/testMeshResolution/2D/mesh1cm', '/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/testMeshResolution/2D/mesh2cm']
    # simPath = ['stopConcAdapTmstp/scat_6-sameDomain/lowCont_lowPe_seed100', 'stopConcAdapTmstp/scat_6-sameDomain/lowCont_seed100', 'stopConcAdapTmstp/scat_6-sameDomain/lowCont_highPe_seed100']
    # simPath = ['stopConcAdapTmstp/scat_6-sameDomain/highCont_lowPe_seed100', 'stopConcAdapTmstp/scat_6-sameDomain/highCont_seed100', 'stopConcAdapTmstp/scat_6-sameDomain/highCont_highPe_seed100']
    # simPath = ['variableMecDisp/varMecDisp3D/lowCont_seed100', 'variableMecDisp/varMecDisp3D/highCont_seed100']
    # simPath = Path('scat_6-sameDomain/lowCont_seed100')    
    # simPath = ['stopConcAdapTmstp/scat_6-sameDomain/lowCont_seed100', 'stopConcAdapTmstp/scat_6-sameDomain/highCont_seed100']
    # simPath = ['stopConcAdapTmstp/scat_5-lowContrast/TS3', 'stopConcAdapTmstp/scat_3-highContrast/TS3']
    # latexFolderPath = Path('/home/pmxep5/OneDrive/Nottingham/Write/Articles/PGSFoam/')
    # saveFolderPath = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/', simPath[i]))
    # homeFolderPath = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/', simPath[i]))

    # simPath = Path('stopConcAdapTmstp/scat_3-highContrast/TS%d' % (FS+i))
    # simPath = Path('stopConcAdapTmstp/scat_5-lowContrast/TS%d' % (FS+i))
    # simPath = Path('stopConcAdapTmstp/scat_7-stochReal/lowCont/TS%d' % (FS+i))
    # latexFolderPath = Path('/home/pmxep5/OneDrive/Nottingham/Write/Articles/PGSFoam/')
    # saveFolderPath = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/', simPath))
    # homeFolderPath = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/', simPath))
    
    simPath = ['realismDegree/Herten7_Stochastic', 'realismDegree/Herten8_Ephesia', 'realismDegree/Herten9_Comunian']
    saveFolderPath = Path(os.path.join('/Users/pmxep5/Git/Overleaf/Thesis/images', simPath[i]))
    homeFolderPath = Path(os.path.join('/Users/pmxep5/Git/Hub/OpenFOAM/PGSFlowTransport/tutorials/Herten/', simPath[i]))                                           
# Parse #######################################################################
    # bashParseLog(sim, FS, homeFolderPath) # Import bashParse.py to use bashParseLog
# parseLog function parses the log file from OpenFOAM and stores the relevant data in different lists
    parseLog(homeFolderPath, s, cl, dd, mvel, c, t, m)
    Kave = []
    if not kval:
        kave.append(1e-12) # In case of no setRandomField the average permeability field value is equal to the value assigned to the uniform permeability field in K.orig
    else:
        for n1, n2 in zip(kval[i], kvol[i]): Kave.append(n1*n2)
        kave.append(sum(Kave)) # Average volumetric permeability [m2]
# parseConstants function parses the transportProperties OpenFOAM input file and stores the hydraulic dispersion coefficient value
    D, Mu, rho, g = parseConstants(homeFolderPath)
# parseFieldsDict function parses the setFieldsDict dictionary and outputs the X distance at which the center of the injection volume is placed (not needed if the injection is on the inlet boundary) 
    topoPath = os.path.join(homeFolderPath, "system/topoSetDict")
    if os.path.isfile(topoPath):
        Xbox = parseSetFieldsDict(homeFolderPath)
    else:
        Xbox = 0
    deltaPx = parseInitialConditions(homeFolderPath)
# Processing ##################################################################
    # Compute statistical parameters, Cumulative Inverse Gaussian and its derivatives    
    mu1, mu1NoUnit, mu2, lam, lamNoUnit, T = processConc(homeFolderPath, dd[i], mvel[i], c[i], t[i], Xbox, s, D, Y, dCnorm, dc, tt, Tadv, tLS, dcLS, dcLSnorm, n)
    y[i] = invGaussianCDF(tt[i], mu1NoUnit, lamNoUnit)
    dY = [(y[i][j+s]-y[i][j])/(tt[i][j+s]-tt[i][j]) for j, val in enumerate(y[i][:-s])] # Smooth derivative
    dYnorm = np.array(dY)/np.array(sum(dY)) # Normalisation of the derivative
    # Store the results in for i-loop into sublists
    dYYnorm[i] = dYnorm
    # Mean advective velocity and Macro Dispersion estimation from statistical moments (Yu 1999)
    v = dd[i][0]/mu1
    MD = [dd[i][0]**2/(2*lam)]
    HD.append(D[0]+MD[0])
    alphaE = MD/v
    alpha = MD[0]/mvel[i][0]
    # van Genuchten approximation of the analytical solution for the ADE
    # yVG.append(invGauVG(np.array(1/mu1NoUnit), np.array(1), np.array(1/(2*lamNoUnit)), np.array(tt[i])))
    # yVG.append(invGauVG(np.array(v), np.array(dd[i][0]), np.array(HD[i]), np.array(t)))
    # yVG.append(invGauVG(np.array(dd[i][0]/mu1NoUnit), np.array(dd[i][0]), np.array(dd[i][0]**2/(2*lamNoUnit)), tt))
    yVG.append(invGauVG(np.array(dd[i][0]/mu1NoUnit), np.array(dd[i][0]), np.array(dd[i][0]**2/(2*lamNoUnit)), np.array(tt[i])))
    # yVG.append(invGauVG(np.array(mvel[i][0]), np.array(dd[i][0]), np.array(HD[i]), np.array(t)))
    # Péclet
    macroPeX = dd[i][0]*mvel[i][0]/D[0]
    macroPeY = dd[i][1]*mvel[i][1]/D[0]
    macroPeZ = dd[i][2]*mvel[i][2]/D[0]
    ndHD.append(HD[i]*Tadv[i]/kave[i])    
    if not cl or cl[i-1][0] == 0:
        cl.append([0, 0, 0]) # In case the permeability field is generated out of OpenFOAM and it has unknown correlation length     
    medPeX = cl[i][0]*mvel[i][0]/D[0]
    medPeY = cl[i][1]*mvel[i][1]/D[0]
    medPeZ = cl[i][2]*mvel[i][2]/D[0]
    # Concentration vector norm
    normC = np.linalg.norm(c[i])
    normY = np.linalg.norm(y[i])
    # Effective permeatbility
    kxeff.append(-mvel[i][0]*Mu[0]/(deltaPx/dd[i][0]-rho[0]*g[2]))
    # Effective permeability / Average volumetric permeability
    kr.append(kxeff[i]/kave[i])
    
# PRINT SECTION ###############################################################
    # From log    
    with open("log") as stats: # It re-opens the log file to print the field statistics 
        for line in stats:
            if "Gaussian Random Field generated" in line:
                print("============= SIMULATION %d =============\n\nGEOSTATISTICAL METRICS" % (i+FS))
                while set(line.split()).isdisjoint(set(["End"])):
                    if "Statistics log" in line:
                        print(line, end = ''),
                        line = stats.readline()
                    print(line, end = '')
                    line = stats.readline()  
    print("Effective longitudinal permeability Kx,eff [m2] = %.15f \n" % kxeff[i])
    print("FLOW METRICS")
    print("Mean Vx Vy Vz = %.9f %.9f %.9f" % (mvel[i][0], mvel[i][1], mvel[i][2]))
    print("Estimated mean Vx = %.9f \n" % (v))  
    print("Péclet: \n  Macro = (%.2f %.2f %.2f) \n  Meso = (%.2f %.2f %.2f) \n" % (macroPeX, macroPeY, macroPeZ, medPeX, medPeY, medPeZ))
    print("TRANSPORT METRICS")
    print("mu_1 = %f \nmu_2 = %f \nInv Gau lambda = %f \nalpha_estimatedMeanVx = %f \nalpha_meanVx = %f \n" % (mu1, mu2, lam, alphaE, alpha))                 
    print("Dispersion: \n  Molecular = %.9f \n  Estimated mechanical = %f \n  Estimated hydrodynamic = %f \n" % (D[0], MD[0], HD[i]))
    print("Computed ||c|| = %f \nEstimated ||c|| = %f \n" % (normC, normY))

# PLOT SECTION ################################################################
font = {'size': 30}
plt.rc('font', **font)

# plt.figure(figsize=(14, 9))
# plt.plot([itm[0] for itm in cl], kr, '*', color='g', label="kx_eff/k_ave")
# plt.plot([itm[0] for itm in cl], ndHD, '*', color='b', label="ndHD")
# plt.xlabel("Lx [m]")
# plt.ylabel("Permeability ratio and non-dimensional hydraulic dispersion [-]")
# plt.legend()
# os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# # plt.savefig(os.path.join(saveFolderPath, "../images/KaveKeff.png"))
# # plt.show()

# plt.figure(figsize=(14, 9))
# plt.plot([itm[0] for itm in cl], kxeff, color='g', label="kx_eff")
# plt.xlabel("Lx [m]")
# plt.ylabel("Effective permeability [m2]")
# plt.legend()
# os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# # plt.savefig(os.path.join(saveFolderPath, "../images/Keff.png"))
# # plt.show()

# plt.figure(figsize=(14, 9))
# plt.plot([itm[0] for itm in cl], HD, color='b', label="hydr disp")
# plt.xlabel("Lx [m]")
# plt.ylabel("Hydraulic dispersion [m2/s]")
# plt.legend()
# os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# # plt.savefig(os.path.join(saveFolderPath, "../images/HD.png"))
# # plt.show()

# plt.figure(figsize=(14, 9))
# plt.plot([itm[0] for itm in cl], HD, color='b', label="hydr disp [m2/s]")
# plt.plot([itm[0] for itm in cl], [itm[0] for itm in mvel], color='r', label="mean Vx  [m/s]")
# plt.xlabel("Lx [m]")
# plt.ylabel("Hydraulic dispersion and mean Vx")
# plt.legend()
# os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# # plt.savefig(os.path.join(saveFolderPath, "../images/HDvsMvel.png"))
# # plt.show()

plt.figure(figsize=(14, 9))
lin = ['-', '-', '-', '-','-', '-', '-', '-', '-', '-']
lab = ['mesh05cm', 'mesh1cm', 'mesh2cm']
# lab = ['mesh1mm', 'mesh 0.5cm', 'mesh 1cm', 'mesh 2cm']
# lab = ['TS1', 'TS2', 'TS3', 'TS4', 'TS5', 'TS6', 'TS7', 'TS8', 'TS9', 'TS10']
# lab = ['Lx = 0.4', 'Lx = 0.6', 'Lx = 0.8', 'Lx = 1.0']
# lab = ['Low Péclet', 'Medium Péclet', 'High Péclet', 'VarMecDisp']
# lab = ['Low k contrast', 'High k contrast']
# lab = ['Dmec = constant', 'Dmec = alpha*V']
col = ['blue', 'orange', 'green', 'red']
# col = ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0']
# col = ['0.0', '0.3', '0.6', '0.9']
# col = ['0.3', '0.6', '0.9']
# col = ['0.15', '0.75']
# col = ['green', 'red']
for j in range(0, sim, interval):
    plt.plot(tLS[j], dcLS[j], ls="%s" % lin[j], color="%s" % col[j], lw=4, label="%s" % lab[j])
    plt.axis([0, 2, 0, max(max(dcLS, key=max))+0.05*max(max(dcLS, key=max))])
    plt.xlabel("T [-]")
    plt.ylabel("$d\overline{c}/dT$ [-]")
    plt.legend()
for j in range(0, sim, interval):
    zoom = plt.axes([.45, .3, .4, .4])
    zoom.plot(tLS[j], dcLS[j], ls="%s" % lin[j], color="%s" % col[j], lw=2)
    zoom.axis([-0.02*max(max(tLS, key=max)), max(max(tLS, key=max)), 0, max(max(dcLS, key=max))+0.05*max(max(dcLS, key=max))])
os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# plt.savefig(os.path.join(saveFolderPath, "../images/PDF.png"))
# plt.savefig(os.path.join(saveFolderPath, "../images/constVarMecDisp.png"))
# plt.show()

plt.figure(figsize=(14, 9))
for j in range(0, sim, interval):
    plt.loglog(tLS[j], dcLS[j], ls="%s" % lin[j], color="%s" % col[j], lw=4, label="%s" % lab[j])
    plt.legend(loc="best", fontsize=30)
    plt.axis([0.1, max(max(tLS, key=max)), 1e-3, max(max(dcLS, key=max))+0.5*max(max(dcLS, key=max))])
    plt.xlabel("T [-]")
    plt.ylabel("$d\overline{c}/dT$ [-]")
plt.tight_layout()
# plt.savefig(os.path.join(latexFolderPath, "images/lowHighCpdf.png"))
# plt.savefig(os.path.join(latexFolderPath, "images/varPeLowC.png"))
# plt.savefig(os.path.join(latexFolderPath, "images/varPeHighC.png"))
# plt.savefig(os.path.join(saveFolderPath, "../images/increasingLxLC.png"))
# plt.savefig(os.path.join(latexFolderPath, "images/increasingLxLC.png"))
# plt.show()

dcLSave = []
tLSave = []
plt.figure(figsize=(14, 9))
for j in range(0, sim, interval):
    plt.loglog(tLS[j], dcLS[j], ls="%s" % lin[j], color="gray", lw=0.5)
    plt.axis([0.1, max(max(tLS, key=max)), 1e-3, max(max(dcLS, key=max))+0.5*max(max(dcLS, key=max))])
    plt.xlabel("T [-]")
    plt.ylabel("$d\overline{c}/dT$ [-]")    
for r in range(0, len(dcLS[0])):
    dcLSave.append(mean([row[r] for row in dcLS]))
    tLSave.append(mean([row[r] for row in tLS]))
plt.loglog(tLSave, dcLSave, color="black", lw=4)
plt.tight_layout()
os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# plt.savefig(os.path.join(latexFolderPath, "images/realisVarLC.png"))
# plt.savefig(os.path.join(saveFolderPath, "../images/increasingLx.png"))
# plt.savefig(os.path.join(saveFolderPath, "../images/logConstVarMecDisp.png"))
# plt.show()

plt.figure(figsize=(14, 9))
for j in range(0, sim, interval):
    # cBoolean = np.logical_and(np.array(dc[j])>1e-3, np.array(dc[j])<1)
    cBoolean = np.logical_and(np.array(c[j])>1e-3, np.array(c[j])<1)
    tThrs = [val for z, val in enumerate(tt[j][:-s]) if cBoolean[z]] # it selects the time only if cBoolean is True
    cThrs = [val for z, val in enumerate(c[j][:-s]) if cBoolean[z]]
    plt.loglog(tThrs, cThrs, ls="%s" % lin[j], color="%s" % col[j], lw=4, label="%s" % lab[j])
    plt.legend(loc="best", fontsize=30)
    plt.xlabel("T [-]")
    plt.ylabel("$\overline{c} [-]$")
    plt.grid(True, which="both")
plt.tight_layout()
# plt.savefig(os.path.join(latexFolderPath, "images/lowHighCcdf.png"))

# plt.figure(figsize=(14, 9))
# plt.plot(tt[i], yVG[i], color='r', label="vanGenuchten")
# plt.plot(tt[i], c[i], color='b', label="sim_BTC")
# plt.xlabel("t [s]")
# plt.ylabel("c [-]")
# plt.legend()
# os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# # plt.savefig(os.path.join(saveFolderPath, "../images/vanGenuchten.png"))
# # plt.show()

# plt.figure(figsize=(14, 9))
# lin = ['-', '-', '-', '-']
# lab = ['Lx = 0.4', 'Lx = 0.6', 'Lx = 0.8', 'Lx = 1.0']
# col = ['green', 'yellow', 'blue', 'orange']
# for j in range(0, sim, interval):
#     plt.loglog(tLS[j], dCnorm[j], ls="%s" % lin[j], color="%s" % col[j], lw=4, label="%s" % lab[j])
#     plt.axis([0.1, 20, 1e-5, 1.6e-3])
#     plt.xlabel("t/t* [-]")
#     plt.ylabel("c [-]")
#     plt.legend()
# os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# plt.savefig(os.path.join(saveFolderPath, "../images/increasingLx.png"))
# # plt.show()

# Draw blank canvas, grids and legend #########################################
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(15, 14))
# font = {'size': 20}
# plt.rc('font', **font)
#plt.tight_layout()

ax[0][0].set_ylabel('$c* [-]$')
ax[1][0].set_ylabel('$dc*/dt* [-]$')
ax[2][0].set_xlabel('$t* [-]$')
ax[2][0].set_ylabel('$1-c* [-]$')
ax[2][1].set_xlabel('$t* [-]$')
ax[2][2].set_xlabel('$t* [-]$')
    
# Plot section ################################################################
for i in [ax[0][0], ax[0][1], ax[0][2], ax[2][0], ax[2][1], ax[2][2]]:
    i.set_ylim([min(min(Y)), max(max(Y))+0.05*max(max(Y))]) # Conc interval (Y-axis)
for i in [ax[1][0], ax[1][1], ax[1][2]]:
    i.set_ylim([1e-5, max(max(dCnorm, key=max))+0.05*max(max(dCnorm, key=max))])
for i in [ax[0][0], ax[0][1]]:    
    i.set_xlim([min(min(tt, key=min)), max(max(tt, key=max))]) # Time interval (X-axis)
for i in [ax[1][0], ax[1][1]]:    
    i.set_xlim([min(min(tt, key=min)), max(max(tt, key=max))]) # Time interval (X-axis)
for i in [ax[2][0], ax[2][1]]:    
    i.set_xlim([min(min(tt, key=min)), max(max(tt, key=max))]) # Time interval (X-axis). Otherwise: min(tt[FS-1:sim])
for i in [ax[0][2], ax[1][2], ax[2][2]]:
    i.set_xlim([1e-1, max(max(tt, key=max))]) # Time interval (X-axis)

for i in range(0, sim, interval):
    for j in range(3): 
        for k in range(3):
            ax[j-1][k-1].grid('on') 
    axs = plt.gca()
    color=next(axs._get_lines.prop_cycler)['color']
    #plt.ticklabel_format(scilimits=[-2, 2])

    # Plot continuous line legend for sim-1 simulations and full legend for the last simulation
    if i+1 < sim:
        mylabelExp = 'TS%d' % (FS+i)
        mylabelIG = '_nolegend_'
        mylabelVG = '_nolegend_'
    elif i+1 == sim:
        mylabelExp = 'TS%d experimental' % (FS+i)
        mylabelIG = 'TS%d cum inv Gau' % (FS+i)
        mylabelVG = 'TS%d van Genuchten' % (FS+i)

    ax[0][0].plot(tt[i], c[i], color=color, label=mylabelExp)
    ax[0][1].semilogy(tt[i], c[i], color=color, label=mylabelExp)
    ax[0][2].loglog(tt[i], c[i], color=color, label=mylabelExp)
    ax[0][0].plot(tt[i], y[i], '--', color=color, label=mylabelIG)
    ax[0][1].semilogy(tt[i], [i for i in y[i]], '--', color=color, label=mylabelIG)
    ax[0][2].loglog(tt[i], [i for i in y[i]], '--', color=color, label=mylabelIG)
    # ax[0][0].plot(tt[i], yVG[i], '-.', color=color, label=mylabelVG)
    # ax[0][1].semilogy(tt[i], yVG[i], '-.', color=color, label=mylabelVG)
    # ax[0][2].loglog(tt[i], yVG[i], '-.', color=color, label=mylabelVG)
    
    ax[1][0].plot(tt[i][:-s], dCnorm[i], color=color, label=mylabelExp)
    ax[1][1].semilogy(tt[i][:-s], dCnorm[i], color=color, label=mylabelExp)
    ax[1][2].loglog(tt[i][:-s], dCnorm[i], color=color, label=mylabelExp)  
    ax[1][0].plot(tt[i][:-s], dYYnorm[i], '--', color=color, label=mylabelIG)
    ax[1][1].semilogy(tt[i][:-s], dYYnorm[i], '--', color=color, label=mylabelIG)
    ax[1][2].loglog(tt[i][:-s], dYYnorm[i], '--', color=color, label=mylabelIG)
    
    ax[2][0].plot(tt[i], [1-i for i in c[i]], color=color, label=mylabelExp)
    ax[2][1].semilogy(tt[i], [1-i for i in c[i]], color=color, label=mylabelExp)
    ax[2][2].loglog(tt[i], [1-i for i in c[i]], color=color, label=mylabelExp)
    ax[2][0].plot(tt[i], 1-y[i], '--', color=color, label=mylabelIG)
    ax[2][1].semilogy(tt[i], [1-i for i in y[i]], '--', color=color, label=mylabelIG)
    ax[2][2].loglog(tt[i], [1-i for i in y[i]], '--', color=color, label=mylabelIG)
    ax[2][0].plot(tt[i], [1-i for i in yVG[i]], '-.', color=color, label=mylabelVG)
    ax[2][1].semilogy(tt[i], [1-i for i in yVG[i]], '-.', color=color, label=mylabelVG)
    ax[2][2].loglog(tt[i], [1-i for i in yVG[i]], '-.', color=color, label=mylabelVG)

ax[0][1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.4), ncol=3)
#ax[0][1].legend(loc='lower center', bbox_to_anchor=(0.5, 0.1), ncol=sim)

# fig.legend(bbox_to_anchor=(0.3,1.1), ncol=4)
# ax.legend(loc=0)
# Adjust visualization, save and show plots
os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# fig.savefig(os.path.join(saveFolderPath, "../images/BTC9axs.png"))
# plt.show()