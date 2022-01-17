#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Process function for OpenFOAM outputs post-processing
"""
import numpy as np
import pandas as pd
#import os

def processConc(path, dd, mvel, c, t, Xbox, s, D, Y, dCnorm, dC, tt, Tadv, dc):
    Ymin = 0 # 1-max(c) # Y minimum plotted value
    # Select significant concentration 
    cBoolean = np.logical_and(np.array(c)>Ymin, np.array(c)<1) 
    Y.append([val for j, val in enumerate(c) if cBoolean[j]])
    # Non-dimensionlise the time
    tadv = (dd[0]-Xbox)/mvel[0] # If punctual injection Xbox needs to be rested from dd 
    Tadv.append(np.array(tadv))
    T = [val/tadv for j, val in enumerate(t) if cBoolean[j]] # it selects the time only if cBoolean is True
    ndT = [val/tadv for j, val in enumerate(t)]
    # tt.append(np.array(ndT))
    
    # CONCENTRATION TIME DERIVATIVE OPTIONS
    n = 1 # Derivative smoothing factor
    sdc = np.array([(c[j+n]-c[j])/(ndT[j+n]-ndT[j]) for j, val in enumerate(c[:-n])]) # Smooth derivative
    sdc = (pd.Series(sdc).rolling(window=n).mean().iloc[s-1:].values)
    dc.append(sdc)
    
    tThrs = [round(val, 11) for z, val in enumerate(ndT[:-s]) if cBoolean[z]] # it selects the time only if cBoolean is True
    dcThrs = [val for z, val in enumerate(sdc) if cBoolean[z]]
    logSpacing = np.logspace(np.log10(min(tThrs)), np.log10(max(tThrs)), 100, endpoint=True)
    logSpacing = np.insert(logSpacing, 0, 0) # it adds 0 at the 0th position
    tLog = []
    dcLog = []
    for k in range(0, len(logSpacing)-1, 1):
        Tlog = []
        dClog = []
        tBoolean = np.logical_and(tThrs>logSpacing[k], tThrs<=logSpacing[k+1])
        Tlog.append([val for j, val in enumerate(tThrs) if tBoolean[j]])
        dClog.append([val for j, val in enumerate(dcThrs) if tBoolean[j]])
        tLog.append(np.mean(Tlog))
        dcLog.append(np.mean(dClog))                                     
    tt.append(np.array(tLog))
    dC.append(dcLog)
    dCnorm.append(dcLog/sum(dcLog))
    # dC.append(dc)
    # dCnorm.append(dc/sum(dc))
    
    # Inverse Gaussian parameters estimation
    c_DT = [] # Needed to perform the dot product
    cDT_t = []
    deltaT = np.array(t[1:])-np.array(t[:-1])
    for n1, n2 in zip(c[:-1], deltaT): c_DT.append(n1*n2) # dot products Matlab equivalent
    for n1, n2 in zip(c_DT, t[:-1]): cDT_t.append(n1*n2) 
    # mu1 is the first order moment of the concentration distribution computed using the integral by parts since 
    # we only know the cdf of the concentration which is the BTC while the pdf is unknown.
    mu1 = c[-1]*t[-1]-np.sum(c_DT)
    mu2 = (c[-1]*(t[-1])**2)-2*np.sum(cDT_t)
    lam = mu1**3/(mu2-mu1**2)
    c_Dt = [] # Needed to perform the dot product
    cDt_t = []
    deltat = np.array(ndT[1:])-np.array(ndT[:-1])
    for n1, n2 in zip(c[:-1], deltat): c_Dt.append(n1*n2) # dot products Matlab equivalent
    for n1, n2 in zip(c_Dt, ndT[:-1]): cDt_t.append(n1*n2)
    mu1NoUnit = c[-1]*ndT[-1]-np.sum(c_Dt)
    mu2NoUnit = (c[-1]*(ndT[-1])**2)-2*np.sum(cDt_t)    
    lamNoUnit = mu1NoUnit**3/(mu2NoUnit-mu1NoUnit**2)
    return mu1, mu1NoUnit, mu2, lam, lamNoUnit, T