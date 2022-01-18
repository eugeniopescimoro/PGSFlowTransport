#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Log-smoothing of BTC noise by log-binning and averaging
"""
import numpy as np

def logAve(ndT, cBoolean, sdc, s, tLS, dcLS, dcLSnorm):
    tLog = [] # List needed to store the new average times which comes after log-binning
    dcLog = [] # List needed to store the new average concentrations which comes after log-binning
    tThrs = [round(val, 10) for z, val in enumerate(ndT[:-s]) if cBoolean[z]] # it selects the time only if cBoolean is True
    dcThrs = [round(val, 10) for z, val in enumerate(sdc) if cBoolean[z]]
    logSpacing = np.logspace(np.log10(min(tThrs)), np.log10(max(tThrs)), 150, endpoint=True) # Log-spaced list whose values range between time boundaries 
    logSpacing = np.insert(logSpacing, 0, 0) # it adds 0 at the 0th position
    for k in range(0, len(logSpacing)-1, 1):
        Tlog = [] # Support list: it collects the times within the Kith log-interval so that an average can be performed
        dClog = [] # Support list: it collects the concentrations within the Kith log-interval so that an average can be performed
        tBoolean = np.logical_and(tThrs>logSpacing[k], tThrs<=logSpacing[k+1])
        Tlog.append([val for j, val in enumerate(tThrs) if tBoolean[j]])
        dClog.append([val for j, val in enumerate(dcThrs) if tBoolean[j]])
        tLog.append(np.mean(Tlog))
        dcLog.append(np.mean(dClog))                                     
    tLS.append(np.array(tLog))
    dcLS.append(dcLog)
    dcLSnorm.append(dcLog/sum(dcLog))