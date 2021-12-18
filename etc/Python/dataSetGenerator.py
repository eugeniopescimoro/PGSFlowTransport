#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 10:25:55 2021
@author: Eugenio Pescimoro, eugenio.pescimoro@gmail.com
Generate custom concentration and time data following cumulative Inverse Gaussian distribution with specified data
"""
import numpy as np
import math
from scipy.stats import norm
import matplotlib.pyplot as plt
from pathlib import Path
from lmfit import Minimizer, Parameters, report_fit
# import os

ti = 1e-10
tf = 80
data = 100000
l = 1
mu = 3
m = 0.1
resultsFolderPath = Path('/data/PGSFlowTransport/tutorials/RESULTS')

t = np.linspace(ti, tf, data)
c = norm.cdf(np.sqrt(l/t)*(t/mu-1))+math.exp(2*l/mu)*norm.cdf(-np.sqrt(l/t)*(t/mu+1))
# os.makedirs(os.path.join(resultsFolderPath, "customInvGauData"), exist_ok = True)
# os.chdir(os.path.join(resultsFolderPath, "customInvGauData"))
# with open('logFake', 'w') as log:
#     for i in range(len(t)):
#         log.write("Adaptive time = %.3f \n" % t[i])    
#         log.write("Flux out = %.3f \n" % c[i])

###############################################################################
# yInvGau = []
# for i in range(len(t)):
#     yInvGau.append(np.sqrt(l/(2*math.pi*t[i]**3))*math.exp(-(l*(t[i]-mu)**2)/(2*mu**2*t[i])))
yInvGau = m*np.sqrt(l/(2*math.pi*t*t*t))*np.exp(-(l*(t-mu)*(t-mu))/(2*mu*mu*t))

###############################################################################
def err_cumInvGau(paramsIG, t, c):
    L = paramsIG['L']
    MU = paramsIG['MU']
    yCumInvGau = norm.cdf(np.sqrt(L/t)*(t/MU-1))+math.exp(2*L/MU)*norm.cdf(-np.sqrt(L/t)*(t/MU+1))
    return np.array(yCumInvGau)-np.array(c)

paramsIG = Parameters()
paramsIG.add('L', value=5, min=0)
paramsIG.add('MU', value=2, min=0)
# do fit, here with the default leastsq algorithm
minner = Minimizer(err_cumInvGau, paramsIG, fcn_args=(t, c))
resultIG = minner.minimize()
# calculate final result
finalIG = c + resultIG.residual
# write error report
report_fit(resultIG)

###############################################################################
def err_parInvGau(paramsParIG, t, dc):
    l = paramsParIG['l']
    mu = paramsParIG['mu']
    m = paramsParIG['m']
    yParInvGau = m*np.sqrt(l/(2*math.pi*t*t*t))*np.exp(-(l*(t-mu)*(t-mu))/(2*mu*mu*t))
    return yParInvGau-dc

r = 1000
paramsParIG = Parameters()
paramsParIG.add('l', value=5, min=0)
paramsParIG.add('mu', value=2, min=0)    
paramsParIG.add('m', value=0.5, min=0)
# do fit, here with the default leastsq algorithm
minnerParIG = Minimizer(err_parInvGau, paramsParIG, fcn_args=(t[0:r], yInvGau[0:r]))
resultParIG = minnerParIG.minimize()    
# calculate final result
finalParIG = yInvGau[0:r] + resultParIG.residual    
# write error report
report_fit(resultParIG)

###############################################################################
c_DT = [] # Needed to perform the dot product
cDT_t = []
deltaT = t[1:]-t[:-1]
for n1, n2 in zip(c[:-1], deltaT): c_DT.append(n1*n2) # dot products Matlab equivalent
for n1, n2 in zip(c_DT, t[:-1]): cDT_t.append(n1*n2) 
mu1 = c[-1]*t[-1]-np.sum(c_DT)
mu2 = (c[-1]*(t[-1])**2)-2*np.sum(cDT_t)
lam = mu1**3/(mu2-mu1**2)
yCIGMom = norm.cdf(np.sqrt(lam/t)*(t/mu1-1))+math.exp(2*l/mu1)*norm.cdf(-np.sqrt(lam/t)*(t/mu1+1))

###############################################################################
# PLOT SECTION
font = {'size': 24}
plt.rc('font', **font)
plt.figure(figsize=(14, 9))
plt.title('Inverse Gaussian estimation tests')
plt.plot(t, yInvGau, color='g', label='Inv Gau: mu=%.1f, lambda=%.1f, m=%.1f' % (mu, l, m))
plt.plot(t, c, color='g', label='Cum IG: mu=%.1f, lambda=%.1f' % (mu, l))
plt.plot(t, finalIG, '--', color='r', label='LSQ: mu=%.1f, lambda=%.1f, m=%.1f' % (minnerParIG.values['mu'], minnerParIG.values['l'], minnerParIG.values['m']))
plt.plot(t, yCIGMom, '--', color='b', label='Moments: mu=%.1f, lambda=%.1f' % (mu1, lam))
plt.xlabel("Dimensionless time t [-]")
plt.ylabel("Concentration c [-]")
plt.legend()