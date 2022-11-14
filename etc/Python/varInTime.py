#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Plot variance and inertial variance from OF text output
"""
import matplotlib.pyplot as plt
import os
from pathlib import Path
import subprocess

varInTimeRoot = '/Users/pmxep5/Git/Hub/OpenFOAM/PGSFlowTransport/tutorials/Herten/varPeclet/lowPeclet/'
varInTimePath = Path(varInTimeRoot)
os.chdir(varInTimePath)
os.makedirs("./LOGs", exist_ok = True)

#cat log | grep 'Time ' | cut -d ' ' -f3 > LOGs/VarTimeStep
subprocess.run(['/bin/bash', '-c', 'cat varInTime | grep \'Time =\' | cut -d\' \' -f3 > LOGs/Time'])
#cat log | grep 'Var=' | cut -d '=' -f4 > LOGs/Var
subprocess.run(['/bin/bash', '-c', 'cat varInTime | grep \'Var=\' | cut -d\'=\' -f4 > LOGs/Var'])
#cat log | grep 'VarConcX=' | cut -d '=' -f4 | cut -d ' ' -f1 > LOGs/VarConcX
subprocess.run(['/bin/bash', '-c', 'cat varInTime | grep \'VarConcX\' | cut -d\'=\' -f4 | cut -d \' \' -f1 > LOGs/VarConcX'])

t = []
with open(Path(os.path.join(varInTimeRoot, "LOGs/Time")), 'r') as Time:
    for line in Time.readlines():
        t.append(float(line))
del t[1]

var = []
with open(Path(os.path.join(varInTimeRoot, 'LOGs/Var')), 'r') as Var:
    for line in Var.readlines():
        var.append(float(line))
        
ivar = []
with open(Path(os.path.join(varInTimeRoot, 'LOGs/VarConcX')), 'r') as Ivar:
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