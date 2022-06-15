#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: non-linear least-square fitting for Van Genuchten analytical solution of ADE with semi-finite domain in 1D  
"""
###############################################################################
#                       EXAMPLE                                               #
###############################################################################
# from numpy import exp, sin

# def residual(variables, x, data, eps_data):
#     """Model a decaying sine wave and subtract data."""
#     amp = variables[0]
#     phaseshift = variables[1]
#     freq = variables[2]
#     decay = variables[3]

#     model = amp * sin(x*freq + phaseshift) * exp(-x*x*decay)

#     return (data-model) / eps_data

# #############################################################################
# from numpy import linspace, random
# from scipy.optimize import leastsq

# # generate synthetic data with noise
# x = linspace(0, 100)
# eps_data = random.normal(size=x.size, scale=0.2)
# data = 7.5 * sin(x*0.22 + 2.5) * exp(-x*x*0.01) + eps_data

# variables = [10.0, 0.2, 3.0, 0.007]
# out = leastsq(residual, variables, args=(x, data, eps_data))

# #############################################################################
# import matplotlib.pyplot as plt

# LSQmodel = out[0][0] * sin(x*out[0][2] + out[0][1]) * exp(-x*x*out[0][3])

# plt.plot(x, data, '*')
# plt.plot(x, LSQmodel)

###############################################################################
#                       APPLICATION                                           #
###############################################################################
import math
import numpy as np
from pathlib import Path
from Parse import parseLog, parseConstants, parseSetFieldsDict
import os
from Process import processConc
from scipy.optimize import leastsq

s = 10 # Smoothness plotting factor: it reads one value every "s"
homeFolderPath = Path('/data/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_3/TS1')
cl, dd, mass, mvel, conc, time = parseLog(homeFolderPath, s)
D, mu, rho, g = parseConstants(homeFolderPath)
topoPath = os.path.join(homeFolderPath, "system/topoSetDict")
if os.path.isfile(topoPath):
    Xbox = parseSetFieldsDict(homeFolderPath)
else:
    Xbox = 0
mu1, mu1NoUnit, mu2, lam, lamNoUnit, t, dCnorm, Y, T = processConc(homeFolderPath, dd, mvel, conc, time, Xbox, s, D)
cBoolean = np.logical_and(np.array(conc)>0, np.array(conc)<1) 
c = [val for i, val in enumerate(conc) if cBoolean[i]]
t = [val for i, val in enumerate(t) if cBoolean[i]]

def err_invGauVG(parameters, t, c):
    u0 = parameters[0]
    X = 2
    D0 = parameters[1]
    
    yVanG = []
    for i in range(len(t)):
        yVanG.append(0.5*math.erfc((X-u0*t[i])/(2*np.sqrt(D0*t[i])))+0.5*math.exp(u0*X/D0)*math.erfc((X+u0*t[i])/(2*np.sqrt(D0*t[i])))+0.5*(2+u0*X/D0+u0**2*t[i]/D0)*math.exp(u0*X/D0)*math.erfc((X+u0*t[i])/(2*np.sqrt(D0*t[i])))-(u0**2*t[i]/(math.pi*D0))**0.5*math.exp(u0*X/D0-(X+u0*t[i])**2/4*D0*t[i]))
    
    return np.array(c)-np.array(yVanG)

parameters = [7.48813e-06, dd[0][0]**2/(2*lam)]
out = leastsq(err_invGauVG, parameters, args=(np.array(t), np.array(c)))

#############################################################################
import matplotlib.pyplot as plt

X = 2
yVanG = []
for i in range(len(t)):
    yVanG.append(0.5*math.erfc((X-out[0][0]*t[i])/(2*np.sqrt(out[0][1]*t[i])))+0.5*math.exp(out[0][0]*X/out[0][1])*math.erfc((X+out[0][0]*t[i])/(2*np.sqrt(out[0][1]*t[i])))+0.5*(2+out[0][0]*X/out[0][1]+out[0][0]**2*t[i]/out[0][1])*math.exp(out[0][0]*X/out[0][1])*math.erfc((X+out[0][0]*t[i])/(2*np.sqrt(out[0][1]*t[i])))-(out[0][0]**2*t[i]/(math.pi*out[0][1]))**0.5*math.exp(out[0][0]*X/out[0][1]-(X+out[0][0]*t[i])**2/4*out[0][1]*t[i]))    

plt.plot(t, c, '*')
plt.plot(t, yVanG)
