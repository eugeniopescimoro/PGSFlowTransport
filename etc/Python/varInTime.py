#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Plot variance and inertial variance from OF text output
"""
import matplotlib.pyplot as plt

t = []
with open('/data/pmxep5-8/PGSFlowTransport/tutorials/stopConcAdapTmstp/scat_5-lowContrast/TS1/Time', 'r') as Time:
    for line in Time.readlines():
        t.append(float(line))

var = []
with open('/data/pmxep5-8/PGSFlowTransport/tutorials/stopConcAdapTmstp/scat_5-lowContrast/TS1/Var', 'r') as Var:
    for line in Var.readlines():
        var.append(float(line))
        
ivar = []
with open('/data/pmxep5-8/PGSFlowTransport/tutorials/stopConcAdapTmstp/scat_5-lowContrast/TS1/VarConcX', 'r') as Ivar:
    for line in Ivar.readlines():
        ivar.append(float(line))

y = []
m = (ivar[2]-ivar[1])/(t[2]-t[1])
for i in range(len(ivar)):
    y.append(m*(t[i]-t[1])+ivar[1]) 

font = {'size': 50}
plt.rc('font', **font)
plt.figure(figsize=(14, 9))
plt.grid(True, which="both")
plt.loglog(t, var, lw=5, color='0.80', label="variance", marker="s", markersize=10)
plt.loglog(t, ivar, lw=5, color='0.40', label='inertial variance', marker="s", markersize=10)
plt.loglog(t, y, lw=5, color='blue', label='ref slope')
plt.legend(loc="best", fontsize=30)

plt.figure(figsize=(14, 9))
plt.grid(True, which="both")
plt.plot(t, var, lw=5, color='0.80', label="variance", marker="s", markersize=10)
plt.plot(t, ivar, lw=5, color='0.40', label='inertial variance', marker="s", markersize=10)
plt.plot(t, y, lw=5, color='blue', label='ref slope')
plt.legend(loc="best", fontsize=30)