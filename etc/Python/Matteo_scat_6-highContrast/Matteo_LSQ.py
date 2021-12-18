#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Constrained least square parameters estimation for the ADE analytical solution in semi-infinite domain provided by Van Genuchten 
Inspired from: https://lmfit.github.io/lmfit-py/parameters.html
"""
###############################################################################
import math
import numpy as np
from pathlib import Path
import os
from lmfit import Minimizer, Parameters, report_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy import special
import subprocess
import re
###############################################################################
s = 10 # Smoothness plotting factor: it reads one value every "s"
mvel = []
time = []
conc = []
dCnorm = []
finalIG = []

###############################################################################

with open("LOGs/logTime") as logTime:
    for line in logTime.readlines():
        time.append(float(line))
with open("LOGs/logFlux") as logFlux:
    for line in logFlux.readlines():
        conc.append(float(line))
with open("LOGs/logVelx") as logVelx:
    for line in logVelx.readlines():
        mvel.append(float(line))

time = time[0:-1:s]
conc = conc[0:-1:s]

Tadv = 2/mvel[0] # If punctual injection Xbox needs to be rested from dd 
ndT = [val/Tadv for j, val in enumerate(time)]
dC = [(conc[j+s]-conc[j])/(ndT[j+s]-ndT[j]) for j, val in enumerate(conc[:-s])] # Smooth derivative
dCnorm.append(np.array(dC)/np.array(sum(dC)))

###############################################################################
def err_InvGau(paramsIG, t, dc):
    l = paramsIG['l']
    mu = paramsIG['mu']
    m = paramsIG['m']
    yInvGau = m*np.sqrt(l/(2*math.pi*t*t*t))*np.exp(-(l*(t-mu)*(t-mu))/(2*mu*mu*t))
    # yInvGau = m*np.sqrt(l/(2*math.pi*t*t*t))*np.exp(-(l*(t-mu)*(t-mu))/(2*mu*mu*t))/sum(dC)
    # yParInvGau = yParInvGau/sum(yParInvGau)
    return yInvGau[:-s]-dc

# LEAST SQUARES METHOD
paramsIG = Parameters()
paramsIG.add('l', value=1)
paramsIG.add('mu', value=3, max=100)
paramsIG.add('m', value=0.1)    
# do fit, here with the default leastsq algorithm
minnerIG = Minimizer(err_InvGau, paramsIG, fcn_args=(np.array(ndT), np.array(dCnorm)))
resultIG = minnerIG.minimize()    
# calculate final result
finalIG = dCnorm + resultIG.residual  

print(">>> PDF FITTING")
report_fit(resultIG)

# PLOT SECTION ################################################################
font = {'size': 24}
plt.rc('font', **font)
plt.figure(figsize=(14, 9))
plt.semilogy(ndT[:-s], dCnorm[0], lw=5, color='blue')
plt.semilogy(ndT[:-s], finalIG[0], linestyle='dotted', lw=5, color='green')
plt.xlim([0, 7]) 
plt.ylim([1e-5, max(max(dCnorm, key=max))+0.05*max(max(dCnorm, key=max))])
plt.xlabel("Dimensionless time t [-]")
plt.ylabel("Concentration c [-]")
plt.legend()