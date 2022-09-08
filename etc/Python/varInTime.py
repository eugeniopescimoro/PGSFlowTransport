#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Plot variance and inertial variance from OF text output
"""
import matplotlib.pyplot as plt

t = []
with open('/Users/pmxep5/Git/Hub/OpenFOAM/PGSFlowTransport/tutorials/stopConcAdapTmstp/scat_5-lowContrast/TS1/Time', 'r') as Time:
    for line in Time.readlines():
        t.append(float(line))

var = []
with open('/Users/pmxep5/Git/Hub/OpenFOAM/PGSFlowTransport/tutorials/stopConcAdapTmstp/scat_5-lowContrast/TS1/Var', 'r') as Var:
    for line in Var.readlines():
        var.append(float(line))
        
ivar = []
with open('/Users/pmxep5/Git/Hub/OpenFOAM/PGSFlowTransport/tutorials/stopConcAdapTmstp/scat_5-lowContrast/TS1/VarConcX', 'r') as Ivar:
    for line in Ivar.readlines():
        ivar.append(float(line))

y = []
m = (ivar[3]-ivar[2])/(t[3]-t[2])
for i in range(len(ivar)):
    y.append(m*(t[i]-t[2])+ivar[2]) 

plt.figure(figsize=(14, 9))
plt.grid()
plt.loglog(t, var, lw=5, color='0.80', label='Low k contrast', marker="s", markersize=10)
plt.loglog(t, ivar, lw=5, color='0.40', label='Low k contrast', marker="s", markersize=10)
plt.loglog(t, y)