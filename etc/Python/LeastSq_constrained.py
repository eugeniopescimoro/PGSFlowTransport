#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Constrained least square parameters estimation for the ADE analytical solution in semi-infinite domain provided by Van Genuchten 
Inspired from: https://lmfit.github.io/lmfit-py/parameters.html
"""
import numpy as np
from pathlib import Path
# from bashParse import bashParseLog, parseLog, parseConstants, parseInitialConditions, parseSetFieldsDict
from Parse import parseLog, parseConstants, parseInitialConditions, parseSetFieldsDict
import os
from Process import processConc
from lmfit import Minimizer, Parameters, report_fit
from InvGau import invGaussianPDF, invGaussianCDF, invGauVG
import matplotlib.pyplot as plt

# INITIALIZATION ##############################################################
sim = 6 # Number of simulations to analyse 
FS = 1 # Number of the First Simulation to analyse
interval = 1 # Interval between increasing simulations
tt = []
c = []
dc = []
t = []
m = []
dCnorm = []
dYcdfmom = []
dYcdflsq = []
Y = []
yVG = [[]]*sim
finalIG = [[]]*sim
finalCIG = [[]]*sim
finalVG = [[]]*sim
yIG = [[]]*sim
yPeakIG = [[]]*sim
finalLogIG = [[]]*sim
peakIG = [[]]*sim
finalPeakIG = [[]]*sim
s = 10 # Smoothness plotting factor: it reads one value every "s"
dd = []
mvel = []
Tadv = []
cl = []
y = []
mass = []
dY = []
dYlsq = []
tLS = []
dcLS = []
dcLSnorm = []
n = 1 # Derivative smoothing factor
diff1 = []
diff2 = []
diff3 = []
diff4 = []

###############################################################################
def err_invGauVG(paramsVG, t, c):
    u0 = paramsVG['u0']
    X = np.array(2)
    D0 = paramsVG['D0']
    yVanG = invGauVG(u0, X, D0, t)
    return yVanG-c

###############################################################################
def err_InvGau(parIG, t, dc):
    l = parIG['l']
    mu = parIG['mu']
    m = 1
    yInvGau = invGaussianPDF(t, l, mu, m)
    return yInvGau[:-s]-dc

###############################################################################
def err_InvGauPeak(paramsPeakIG, t, dc):
    l = paramsPeakIG['l']
    mu = paramsPeakIG['mu']
    m = paramsPeakIG['m']
    yInvGau = invGaussianPDF(t, l, mu, m)
    return yInvGau[:-s]-dc

###############################################################################
def err_cumInvGau(paramsCIG, t, c):
    l = paramsCIG['l']
    mu = paramsCIG['mu']
    yCumInvGau = invGaussianCDF(t, mu, l)
    return yCumInvGau-c

# ###############################################################################
# def err_LogInvGau(paramsLogIG, Logt, Logdc):
#     l = paramsLogIG['l']
#     mu = paramsLogIG['mu']
#     yLogInvGau = np.log(np.sqrt(l/(2*math.pi*Logt*Logt*Logt))*np.exp(-(l*(Logt-mu)*(Logt-mu))/(2*mu*mu*Logt)))
#     # yInvGau = yInvGau/sum(yInvGau)
#     return yLogInvGau-Logdc

# LOOP THROUGH THE SIMULATIONS ################################################ 
for i in range(0, sim, interval):
# Paths 
    # simPath = ['variableMecDisp/varMecDisp3D/lowCont_seed100', 'variableMecDisp/varMecDisp3D/highCont_seed100']
    # simPath = ['scat_6-sameDomain/lowCont_lowPe_seed100', 'scat_6-sameDomain/lowCont_seed100', 'scat_6-sameDomain/lowCont_highPe_seed100']
    # simPath = ['scat_6-sameDomain/highCont_lowPe_seed100', 'scat_6-sameDomain/highCont_seed100', 'scat_6-sameDomain/highCont_highPe_seed100']
    simPath = ['scat_6-sameDomain/lowCont_lowPe_seed100', 'scat_6-sameDomain/lowCont_seed100', 'scat_6-sameDomain/lowCont_highPe_seed100', 'scat_6-sameDomain/highCont_lowPe_seed100', 'scat_6-sameDomain/highCont_seed100', 'scat_6-sameDomain/highCont_highPe_seed100']
    # simPath = ['stopConcAdapTmstp/scat_6-sameDomain/highCont_lowPe_seed100']
    # simPath = ['scat_3-highContrast/TS3']
    # simPath = ['stopConcAdapTmstp/scat_3-highContrast/TS4']
    # simPath = ['scat_5-lowContrast/TS3']
    # simPath = Path('scat_6-sameDomain/lowCont_seed100')
    # simPath = ['scat_5-lowContrast/TS1', 'scat_5-lowContrast/TS2', 'scat_5-lowContrast/TS3', 'scat_5-lowContrast/TS4', 'scat_3-highContrast/TS1', 'scat_3-highContrast/TS2', 'scat_3-highContrast/TS3', 'scat_3-highContrast/TS4']
    latexFolderPath = Path('/home/pmxep5/OneDrive/Nottingham/Write/Articles/PGSFoam/')
    saveFolderPath = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/', simPath[i]))
    homeFolderPath = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/', simPath[i]))
    # simPath = Path('TS%d' % (FS+i))
    # saveFolderPath = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/variableMecDisp/varMecDisp3D', simPath))
    # homeFolderPath = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/', simPath))
    # homeFolderPath = Path('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_3/TS1')
    # homeFolderPath = Path('/data/PGSFlowTransport/tutorials/RESULTS/scat_5-lowContrast/TS%d' % (FS+i))
    # homeFolderPath = Path('/data/PGSFlowTransport/tutorials/RESULTS/stopConc_5-lowContrast/TS%d' % (FS+i))
    # homeFolderPath = Path('/data/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_3/TS%d' % i)
    # homeFolderPath = Path('/data/PGSFlowTransport/tutorials/RESULTS/Herten/Herten9')
    # homeFolderPath = Path('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_3/TS%d' % (i+FS))

# Parse #######################################################################
    # bashParseLog(sim, FS, homeFolderPath) # TO BE UNCOMMENTED WHEN USING BASHPARSE !!

# parseLog function parses the log file from OpenFOAM and stores the relevant data in different lists
    kvol, kval = parseLog(homeFolderPath, s, cl, dd, mvel, c, t, m)
    
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
    mu1, mu1NoUnit, mu2, lam, lamNoUnit, T = processConc(homeFolderPath, dd[i], mvel[i], c[i], t[i], Xbox, s, D, Y, dCnorm, dc, tt, Tadv, tLS, dcLS, dcLSnorm, n) # Compute statistical parameters, Cumulative Inverse Gaussian and its derivatives
    
    # MOMENTS METHOD
    y.append(invGaussianCDF(tt[i], mu1NoUnit, lamNoUnit)) # CDF estimated with parameters computed on non-dimensional times
    dY.append(np.array([(y[i][j+s]-y[i][j])/(tt[i][j+s]-tt[i][j]) for j, val in enumerate(y[i][:-s])])) # Smooth derivative
    dYcdfmom.append(dY[i]/sum(dY[i])) # Normalisation of the derivative
    v = dd[i][0]/mu1
    MD = dd[i][0]**2/(2*lam)

    # LEAST SQUARES METHOD
    # Van Genuchten
    # paramsVG = Parameters()
    # paramsVG.add('u0', value=1e-06, min=0)
    # # paramsVG.add('X', value=2, min=0)
    # paramsVG.add('D0', value=1e-7, min=0)
    # # do fit, here with the default leastsq algorithm
    # minner = Minimizer(err_invGauVG, paramsVG, fcn_args=(np.array(tt[i]), c[i]))
    # resultVG = minner.minimize()
    # # calculate final result
    # finalVG[i] = c[i] + resultVG.residual

    # Inverse Gaussian
    paramsIG = Parameters()
    paramsIG.add('l', value=lam, min=0)
    paramsIG.add('mu', value=mu1, min=0)
    # paramsIG.add('m', value=1)    
    minnerIG = Minimizer(err_InvGau, paramsIG, fcn_args=(np.array(t[i]), dc[i]))
    resultIG = minnerIG.minimize()
    vIG = dd[i][0]/minnerIG.values['mu']
    MDig = dd[i][0]**2/(2*minnerIG.values['l'])
    paramsIGnoun = Parameters()
    paramsIGnoun.add('l', value=lamNoUnit, min=0)
    paramsIGnoun.add('mu', value=mu1NoUnit, min=0)
    # paramsIGnoun.add('m', value=1)    
    minnerIGnoun = Minimizer(err_InvGau, paramsIGnoun, fcn_args=(np.array(tLS[i]), dcLS[i][:-s]))
    resultIGnoun = minnerIGnoun.minimize()    
    # finalIG[i] = dc[i] + resultIG.residual   
    finalIG[i] = invGaussianPDF(np.array(tt[i][:-s]), minnerIGnoun.values['l'], minnerIGnoun.values['mu'], 1)# minnerIGnoun.values['m'])
    mass.append(np.trapz(dc[i], x=tt[i][:-s])) # Check on total mass -> should be around max(c[i])
    mu1IG = abs(dd[i][0]/vIG - mu1)/mu1*100
    lamIG = abs(dd[i][0]**2/(2*MDig) - lam)/lam*100   
    
    # Cumulative Inverse Gaussian
    paramsCIG = Parameters()
    paramsCIG.add('l', value=lam, min=0)
    paramsCIG.add('mu', value=mu1, min=0)    
    minnerCIG = Minimizer(err_cumInvGau, paramsCIG, fcn_args=(np.array(t[i]), c[i]))
    resultCIG = minnerCIG.minimize()
    vCIG = dd[i][0]/minnerCIG.values['mu']
    MDcig = dd[i][0]**2/(2*minnerCIG.values['l'])
    paramsCIGnoun = Parameters()
    paramsCIGnoun.add('l', value=lamNoUnit, min=0)
    paramsCIGnoun.add('mu', value=mu1NoUnit, min=0) 
    minnerCIGnoun = Minimizer(err_cumInvGau, paramsCIGnoun, fcn_args=(np.array(tt[i]), c[i]))
    resultCIGnoun = minnerCIGnoun.minimize()    
    # finalCIG[i] = c[i] + resultCIG.residual   
    finalCIG[i] = invGaussianCDF(np.array(tt[i]), minnerCIGnoun.values['l'], minnerCIGnoun.values['mu'])
    mu1CIG = abs(dd[i][0]/vCIG - mu1)/mu1*100
    lamCIG = abs(dd[i][0]**2/(2*MDcig) - lam)/lam*100   
 
    # paramsLogIG = Parameters()
    # paramsLogIG.add('l', value=1, min=0)
    # paramsLogIG.add('mu', value=3, min=0)    
    # logDist = [1, 10, 100, 1000, 10000]
    # # do fit, here with the default leastsq algorithm
    # minnerLogIG = Minimizer(err_LogInvGau, paramsLogIG, fcn_args=(np.array([tt[i][j] for j in logDist]), np.array([np.log(dCCnorm[i][j]) for j in logDist])))
    # resultLogIG = minnerLogIG.minimize()    
    # # calculate final result
    # finalLogIG[i] = np.array([dCCnorm[i][j] for j in logDist]) + resultLogIG.residual    
    # # write error report
    # report_fit(resultLogIG)

    # Peak of the Inverse Gaussian
    peakSrt = 50
    peakEnd = 350
    paramsPeakIG = Parameters()
    paramsPeakIG.add('l', value=lam, min=0)
    paramsPeakIG.add('mu', value=mu1, min=0)
    # paramsPeakIG.add('m', value=1, min=0)  
    minnerPeakIG = Minimizer(err_InvGau, paramsPeakIG, fcn_args=(np.array(t[i][peakSrt:peakEnd+s]), dc[i][peakSrt:peakEnd]))
    resultPeakIG = minnerPeakIG.minimize()    
    vPeakIG = dd[i][0]/minnerPeakIG.values['mu']
    MDpeakIG = dd[i][0]**2/(2*minnerPeakIG.values['l'])
    paramsPeakIGnoun = Parameters()
    paramsPeakIGnoun.add('l', value=lamNoUnit, min=0)
    paramsPeakIGnoun.add('mu', value=mu1NoUnit, min=0)
    # paramsPeakIGnoun.add('m', value=1, min=0)    
    minnerPeakIGnoun = Minimizer(err_InvGau, paramsPeakIGnoun, fcn_args=(np.array(tt[i][peakSrt:peakEnd+s]), dc[i][peakSrt:peakEnd]))
    resultPeakIGnoun = minnerPeakIGnoun.minimize()    
    # peakIG[i] = dCnorm[i][peakSrt:peakEnd] + resultPeakIG.residual
    finalPeakIG[i] = invGaussianPDF(np.array(tt[i]), minnerPeakIGnoun.values['l'], minnerPeakIGnoun.values['mu'], 1)# minnerPeakIGnoun.values['m'])
    mu1PeakIG = abs(dd[i][0]/vPeakIG - mu1)/mu1*100
    lamPeakIG = abs(dd[i][0]**2/(2*MDpeakIG) - lam)/lam*100    

    # LSQ DERIVATIVE
    dYlsq.append(np.array([(finalCIG[i][peakEnd+s]-finalCIG[i][peakEnd])/(tt[i][peakEnd+s]-tt[i][peakEnd]) for peakEnd, val in enumerate(finalCIG[i][:-s])])) # Smooth derivative
    dYcdflsq.append(dYlsq[i]/sum(dYlsq[i])) # Normalisation of the derivative
    
    # LSQ INTEGRAL 
    yIG[i] = np.cumsum(finalIG[i])
    yPeakIG[i] = np.cumsum(finalPeakIG[i])

    diff1.append(np.linalg.norm(dc[i]-dY[i])/np.linalg.norm(dc[i])*100)
    diff2.append(np.linalg.norm(dc[i]-finalIG[i])/np.linalg.norm(dc[i])*100)
    diff3.append(np.linalg.norm(dc[i]-dYlsq[i])/np.linalg.norm(dc[i])*100)
    diff4.append(np.linalg.norm(dc[i]-finalPeakIG[i][:-s])/np.linalg.norm(dc[i])*100)

# PRINT SECTION ###############################################################
    print("\n============= SIMULATION %d =============" % (i+FS))
    print(">>> PDF FITTING")
    report_fit(resultIGnoun)
    print(">>> CDF FITTING")
    report_fit(resultCIGnoun)
    print(">>> PEAK PDF FITTING")
    report_fit(resultPeakIGnoun)  # write error report
    # print(">>> VAN GENUCHTEN FITTING")
    # report_fit(resultVG)
    print("\n||Experimental c - Estimated c|| [%%]: \nMoments on CDF %f \nLSQ on PDF %f \nLSQ on CDF %f \nLSQ on PDF Peak %f" % (diff1[i], diff2[i], diff3[i], diff4[i]))
    print("\nEstimated velocity and dispersion: \nMoments on CDF %.9E %.9E \nLSQ on PDF %.9E %.9E \nLSQ on CDF %.9E %.9E \nLSQ on PDF Peak %.9E %.9E" % (v, MD, vIG, MDig, vCIG, MDcig, vPeakIG, MDpeakIG))
    print("\nEstimated velocity discrepancy [%%]: \nMoments on CDF %.9f \nLSQ on PDF %.9f \nLSQ on CDF %.9f \nLSQ on PDF Peak %.9f " % (abs(mvel[0][0]-v)/v*100, abs(mvel[0][0]-vIG)/vIG*100, abs(mvel[0][0]-vCIG)/vCIG*100, abs(mvel[0][0]-vPeakIG)/vPeakIG*100))
    print("\nEstimated dispersion discrepancy [%%]: \nMoments on CDF %.9f \nLSQ on PDF %.9f \nLSQ on CDF %.9f \nLSQ on PDF Peak %.9f " % (abs(mvel[0][0]*cl[0][0]-MD)/MD*100, abs(mvel[0][0]*cl[0][0]-MDig)/MDig*100, abs(mvel[0][0]*cl[0][0]-MDcig)/MDcig*100, abs(mvel[0][0]*cl[0][0]-MDpeakIG)/MDpeakIG*100))
    print("\nEstimated statistical moments: \nLSQ on PDF %.9f %.9f \nLSQ on CDF %.9f %.9f \nLSQ on PDF Peak %.9f %.9f" % (mu1IG, lamIG, mu1CIG, lamCIG, mu1PeakIG, lamPeakIG))
    
# PLOT SECTION ################################################################

# FIGURE 1 ####################################################################
# Draw blank canvas, grids and legend #########################################
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(14, 9))
font = {'size': 50}
plt.rc('font', **font)
#plt.tight_layout()

ax[0][0].set_ylabel('$c* [-]$')
ax[1][0].set_ylabel('$dc*/dt* [-]$')
ax[2][0].set_ylabel('$1-c* [-]$')
ax[2][0].set_xlabel('$t* [-]$')
ax[2][1].set_xlabel('$t* [-]$')
ax[2][2].set_xlabel('$t* [-]$')

for i in [ax[0][0], ax[0][1], ax[0][2], ax[2][0], ax[2][1], ax[2][2]]:
    i.set_ylim([min(min(Y)), max(max(Y))+0.1*max(max(Y))]) # Conc interval (Y-axis)
for i in [ax[1][0], ax[1][1], ax[1][2]]:
    i.set_ylim([1e-5, max(max(dc, key=max))+0.1*max(max(dc, key=max))])
for i in [ax[0][0], ax[0][1]]:    
    i.set_xlim([min(min(tt, key=min)), max(max(tt, key=max))]) # Time interval (X-axis)
for i in [ax[1][0], ax[1][1]]:    
    i.set_xlim([min(min(tt, key=min)), max(max(tt, key=max))]) # Time interval (X-axis)
for i in [ax[2][0], ax[2][1]]:    
    i.set_xlim([min(min(tt, key=min)), max(max(tt, key=max))]) # Time interval (X-axis). Otherwise: min(tt[FS-1:sim])
for i in [ax[0][2], ax[1][2], ax[2][2]]:
    i.set_xlim([1e-1, max(max(tt, key=max))]) # Time interval (X-axis)

for i in range(sim):

    for j in range(3):
        for z in range(3):
            ax[j-1][z-1].grid('on')
    axs = plt.gca()
    color=next(axs._get_lines.prop_cycler)['color']
    #plt.ticklabel_format(scilimits=[-2, 2])

    if i+1 < sim:
        mylabelExp = 'TS%d' % (FS+i)
        mylabelLsq = '_nolegend_'
        mylabelMom = '_nolegend_'
        mylabeldCdT = '_nolegend_'
        mylabelPeak = '_nolegend_'
    elif i+1 == sim:
        mylabelExp = '_nolegend_'
        mylabeldCdT = 'TS%d LSQ PDF' % (FS+i)
        mylabelMom = 'TS%d moments CDF' % (FS+i)
        mylabelLsq = 'TS%d LSQ CDF' % (FS+i)
        mylabelPeak = 'TS%d LSQ Peak PDF' % (FS+i)

    ax[0][0].plot(tt[i], c[i], color=color, label=mylabelExp)
    ax[0][1].semilogy(tt[i], c[i], color=color, label=mylabelExp)
    ax[0][2].loglog(tt[i], c[i], color=color, label=mylabelExp)
    ax[0][0].plot(tt[i], finalCIG[i], '-.', color=color, label=mylabelLsq)
    ax[0][1].semilogy(tt[i], finalCIG[i], '-.', color=color, label=mylabelLsq)
    ax[0][2].loglog(tt[i], finalCIG[i], '-.', color=color, label=mylabelLsq)
    ax[0][0].plot(tt[i], y[i], '--', color=color, label=mylabelMom)
    ax[0][1].semilogy(tt[i], [i for i in y[i]], '--', color=color, label=mylabelMom)
    ax[0][2].loglog(tt[i], [i for i in y[i]], '--', color=color, label=mylabelMom)
    ax[0][0].plot(tt[i][:-s], yIG[i], linestyle='dotted', color=color, label=mylabeldCdT)
    ax[0][1].semilogy(tt[i][:-s], yIG[i], linestyle='dotted', color=color, label=mylabeldCdT)
    ax[0][2].loglog(tt[i][:-s], yIG[i], linestyle='dotted', color=color, label=mylabeldCdT)
    ax[0][0].plot(tt[i], yPeakIG[i], linestyle='dotted', color='red', label=mylabelPeak)
    ax[0][1].semilogy(tt[i], yPeakIG[i], linestyle='dotted', color='red', label=mylabelPeak)
    ax[0][2].loglog(tt[i], yPeakIG[i], linestyle='dotted', color='red', label=mylabelPeak)    
    
    ax[1][0].plot(tt[i][:-s], dc[i], color=color, label=mylabelExp)
    ax[1][1].semilogy(tt[i][:-s], dc[i], color=color, label=mylabelExp)
    ax[1][2].loglog(tt[i][:-s], dc[i], color=color, label=mylabelExp)  
    ax[1][0].plot(tt[i][:-s], dYlsq[i], '-.', color=color, label=mylabelLsq)
    ax[1][1].semilogy(tt[i][:-s], dYlsq[i], '-.', color=color, label=mylabelLsq)
    ax[1][2].loglog(tt[i][:-s], dYlsq[i], '-.', color=color, label=mylabelLsq)
    ax[1][0].plot(tt[i][:-s], dY[i], '--', color=color, label=mylabelMom)
    ax[1][1].semilogy(tt[i][:-s], dY[i], '--', color=color, label=mylabelMom)
    ax[1][2].loglog(tt[i][:-s], dY[i], '--', color=color, label=mylabelMom)
    ax[1][0].plot(tt[i][:-s], finalIG[i], linestyle='dotted', color=color, label=mylabeldCdT)
    ax[1][1].semilogy(tt[i][:-s], finalIG[i], linestyle='dotted', color=color, label=mylabeldCdT)
    ax[1][2].loglog(tt[i][:-s], finalIG[i], linestyle='dotted', color=color, label=mylabeldCdT)
    ax[1][0].plot(tt[i], finalPeakIG[i], linestyle='dotted', color='red', label=mylabeldCdT)
    ax[1][1].semilogy(tt[i], finalPeakIG[i], linestyle='dotted', color='red', label=mylabeldCdT)
    ax[1][2].loglog(tt[i], finalPeakIG[i], linestyle='dotted', color='red', label=mylabeldCdT)
    # ax[1][0].plot(np.array([tt[i][j] for j in logDist]), finalLogIG[i], linestyle='dotted', color='red')
    
    ax[2][0].plot(tt[i], [1-i for i in c[i]], color=color, label=mylabelExp)
    ax[2][1].semilogy(tt[i], [1-i for i in c[i]], color=color, label=mylabelExp)
    ax[2][2].loglog(tt[i], [1-i for i in c[i]], color=color, label=mylabelExp)
    ax[2][0].plot(tt[i], [1-i for i in finalCIG[i]], '-.', color=color, label=mylabelLsq)
    ax[2][1].semilogy(tt[i], [1-i for i in finalCIG[i]], '-.', color=color, label=mylabelLsq)
    ax[2][2].loglog(tt[i], [1-i for i in finalCIG[i]], '-.', color=color, label=mylabelLsq)
    ax[2][0].plot(tt[i], 1-y[i], '--', color=color, label=mylabelMom)
    ax[2][1].semilogy(tt[i], [1-i for i in y[i]], '--', color=color, label=mylabelMom)
    ax[2][2].loglog(tt[i], [1-i for i in y[i]], '--', color=color, label=mylabelMom)
    ax[2][0].plot(tt[i][:-s], [1-i for i in yIG[i]], linestyle='dotted', color=color, label=mylabeldCdT)
    ax[2][1].semilogy(tt[i][:-s], [1-i for i in yIG[i]], linestyle='dotted', color=color, label=mylabeldCdT)
    ax[2][2].loglog(tt[i][:-s], [1-i for i in yIG[i]], linestyle='dotted', color=color, label=mylabeldCdT)

ax[0][1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.35), ncol=4)

os.makedirs(os.path.join(saveFolderPath, "../images"), exist_ok = True)
# fig.savefig(os.path.join(saveFolderPath, "../images/MomVsLsq.png"))
# plt.show()

# FIGURE 2 ####################################################################
# font = {'size': 24}
# plt.rc('font', **font)
plt.figure(figsize=(14, 9))

mylabelExp = 'Simulation data'
mylabelMom = 'Moments on CDF'
mylabeldCdT = 'LSQ on PDF'
mylabelLsq = 'LSQ on CDF'
mylabelPeak = 'LSQ on PDF Peak'

cLSBoolean = np.logical_and(np.array(dcLS[i])>1e-6, np.array(dcLS[i])<1)
# cBoolean = np.logical_and(np.array(dc[i])>1e-4, np.array(dc[i])<1)
plt.loglog(tLS[i][tLS[i]*cLSBoolean != 0], np.array(dcLS[i])[dcLS[i]*cLSBoolean != 0], lw=5, color='0.85', label=mylabelExp)
plt.loglog(tt[i][:-s], dY[i], '--', lw=5, color='0.55', label=mylabelMom)
plt.loglog(tt[i][:-s], finalIG[i], linestyle='dotted', lw=5, color='0.55', label=mylabeldCdT)
plt.loglog(tt[i][:-s], dYlsq[i], '-.', lw=5, color='0.55', label=mylabelLsq)
plt.loglog(tt[i], finalPeakIG[i], linestyle='dotted', lw=5, color='0.35', label=mylabelPeak)
# plt.plot(tt[i], c[i], 'g', label="Experimental data")
# plt.plot(tt[i], finalVG[i], '--', color='r', label="Van Genuchten LSQ")
# plt.plot(tt[i], finalCIG[i], '--', color='b', label="Cum Inv Gau LSQ")

# plt.xlim([min(tt[i][:-s][tt[i][:-s]*cBoolean != 0])-0.1*min(tt[i][:-s][tt[i][:-s]*cBoolean != 0]), max(tt[i][:-s][tt[i][:-s]*cBoolean != 0])]) 
# plt.ylim([min(dc[i][dc[i]*cBoolean != 0]), max(dc[i][dc[i]*cBoolean != 0])+0.05*max(dc[i][dc[i]*cBoolean != 0])])
plt.axis([0.1, max(max(tLS, key=max)), 1e-3, max(max(dcLS, key=max))+0.5*max(max(dcLS, key=max))])

plt.xlabel("T [-]")
plt.ylabel("$d\overline{c}/dT$ [-]")    
plt.legend(fontsize=30)
plt.tight_layout()

os.makedirs(os.path.join(saveFolderPath, "images"), exist_ok = True)
# plt.savefig(os.path.join(latexFolderPath, "images/BTCInterp_lowC_semiLog.png"))
plt.savefig(os.path.join(latexFolderPath, "images/BTCInterp_highC_semiLog.png"))
# plt.savefig(os.path.join(saveFolderPath, "images/BTCInterp_semiLog.png"))
plt.show()

# # FIGURE 3A ###################################################################
# Lx = []
# plt.figure(figsize=(14, 9))
# for i in range(int(sim/2)):
#     Lx.append(cl[0:sim][i][0])
# plt.plot(Lx, diff1[0:4], lw=5, color='0.80', label='Low k contrast')
# plt.plot(Lx, diff1[4:8], lw=5, color='0.4', label='High k contrast')
# plt.xlabel("$\lambda_x [-]$")
# plt.ylabel("e [%]")
# plt.legend(fontsize=30)
# plt.tight_layout()
# plt.savefig(os.path.join(latexFolderPath, "images/LxVsErr.png"))

# FIGURE 3B ###################################################################
lowContPe = [8e2, 8e3, 8e4]
highContPe = [6e3, 6e4, 6e5]
plt.figure(figsize=(14, 9))
plt.semilogx(lowContPe, diff1[0:3], lw=5, color='0.80', label='Low k contrast')
plt.semilogx(highContPe, diff1[3:6], lw=5, color='0.4', label='High k contrast')
plt.xlabel("Pe [-]")
plt.ylabel("e [%]")
plt.legend(fontsize=30)
plt.tight_layout()
plt.savefig(os.path.join(latexFolderPath, "images/errPe.png"))