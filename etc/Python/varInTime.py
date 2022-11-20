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
import numpy as np

# N.B.: varInTime is the output of the postProcess command 
# "mpirun --oversubscribe -np 96 postProcess -dict system/fieldMetricsDict -field c -parallel >> varInTime 2>&1"
# run on the server where concentration fields for different time steps are stored and the following functionObject in the system/fieldMetricsDict
# functions
# {
#     concentrationMetrics
#     {
#         type        fieldMetrics;
#         libs        ("libfieldMetricsFunctionObject.so");
#         operations  ("mean" "meanConc");
#         nRadii      100;
#         maxDist     1.0;
#         fieldName   c;
#     }
# }
# Transfer from the server to local machine with "myscp /data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/Herten/varPeclet/highPeclet/varInTime ."

sim = 3
interval = 1
t = []
var = []
ivar = []

simPath = ['injectionArea/localInjection/', 'injectionArea/wellInjection/', 'injectionArea/wallInjection/']
#simPath = ['varPeclet/lowPeclet/', 'varPeclet/mediumPeclet/', 'varPeclet/highPeclet/']
#simPath = ['varMecDisp/varMecDisp1e-3alpha/', 'varMecDisp/varMecDisp1e-4alpha/', 'varMecDisp/varMecDisp1e-5alpha/']
#simPath = ['realismDegree/Herten7_Stochastic/', 'realismDegree/Herten8_Ephesia/', 'realismDegree/Herten9_Comunian/']
latexFolderPath = Path('/Users/pmxep5/Git/Overleaf/Thesis/')

for i in range(0, sim, interval):
    time = []
    varInTimePath = os.path.join('/Users/pmxep5/Git/Hub/OpenFOAM/PGSFlowTransport/tutorials/Herten/', simPath[i])
    os.chdir(Path(varInTimePath))
    os.makedirs("./LOGs", exist_ok = True)
    #cat log | grep 'Time ' | cut -d ' ' -f3 > LOGs/VarTimeStep
    subprocess.run(['/bin/bash', '-c', 'cat varInTime | grep \'Time =\' | cut -d\' \' -f3 > LOGs/Time'])
    #cat log | grep 'Var=' | cut -d '=' -f4 > LOGs/Var
    subprocess.run(['/bin/bash', '-c', 'cat varInTime | grep \'Var=\' | cut -d\'=\' -f4 > LOGs/Var'])
    #cat log | grep 'VarConcX=' | cut -d '=' -f4 | cut -d ' ' -f1 > LOGs/VarConcX
    subprocess.run(['/bin/bash', '-c', 'cat varInTime | grep \'VarConcX\' | cut -d\'=\' -f4 | cut -d \' \' -f1 > LOGs/VarConcX'])
    with open(Path(os.path.join(varInTimePath, "LOGs/Time")), 'r') as Time:
        for line in Time.readlines():
            time.append(float(line))
    del time[1]
    t.append(np.array(time))
    
    variance = []
    with open(Path(os.path.join(varInTimePath, 'LOGs/Var')), 'r') as Var:
        for line in Var.readlines():
            variance.append(float(line))
    var.append(np.array(variance))
    
    inertialVar = []
    with open(Path(os.path.join(varInTimePath, 'LOGs/VarConcX')), 'r') as Ivar:
        for line in Ivar.readlines():
            inertialVar.append(float(line))
    ivar.append(np.array(inertialVar))

y = []
z = []
a = 1e-6
b = 2
for z in range(len(ivar)):
    y.append(a*t[z]**b)

font = {'size': 50}
plt.rc('font', **font)

lin = ['-', '-.', '--']
col = ['blue', 'orange', 'green', 'red']
#lab = ['Dm=1e-9', 'Dm=1e-10', 'Dm=1e-11']
#lab = ['alpha=1e-3', 'alpha=1e-4', 'alpha=1e-5']
#lab = ['Stochastic', 'Conditiones', 'Realism']
lab = ['Local', 'Well', 'Wall']

plt.figure(figsize=(14, 9))
plt.grid(True, which="both")
for j in range(0, sim, interval):
    plt.loglog(t[j], var[j], ls="%s" % lin[j], lw=5, color="%s" % col[j], label="%s" % lab[j], marker="s", markersize=10)
    # plt.loglog(t, ivar, lw=5, color='0.40', label='inertial variance', marker="s", markersize=10)
    # plt.loglog(t, y, lw=5, color='blue', label='ref slope')
plt.legend(loc="best", fontsize=30)
#plt.savefig(os.path.join(latexFolderPath, "images/PecletViT.png"))
#plt.savefig(os.path.join(latexFolderPath, "images/mecDispViT.png"))
#plt.savefig(os.path.join(latexFolderPath, "images/realDegViT.png"))
plt.savefig(os.path.join(latexFolderPath, "images/injAreaViT.png"))

plt.figure(figsize=(14, 9))
plt.grid(True, which="both")
# N.M.B. PORCO DEMONIO INIZIALIZZARE CON divar = [[]]*sim E' UNA CAGATA PERCHE' CREA TRE SUBLISTS IDENTICHE !!! 
divar = [[] for i in range(sim)] 
for j in range(0, sim, interval):
    # plt.plot(t, var, lw=5, color='0.80', label="variance", marker="s", markersize=10)
    plt.loglog(t[j], ivar[j], ls="%s" % lin[j], lw=5, color="%s" % col[j], label="%s" % lab[j], marker="s", markersize=10)
    for k in range(0, len(ivar[j])-1): # Find the significative range of variances to be plotted (before the plateau)
        divar[j].append((ivar[j][k+1]-ivar[j][k])/(t[j][k+1]-t[j][k]))
plt.loglog(t[j], y[j], ls="-", lw=5, color="black", label='y=a*x^b with a=%f b=%f' % (a, b))
plt.axis([min(min(t, key=min)), 1*10**6, min(min(ivar, key=min)), max(max(ivar, key=max))])
plt.legend(loc="best", fontsize=30)
#plt.savefig(os.path.join(latexFolderPath, "images/PecletIViT.png"))
#plt.savefig(os.path.join(latexFolderPath, "images/mecDispIViT.png"))
#plt.savefig(os.path.join(latexFolderPath, "images/realDegIViT.png"))
plt.savefig(os.path.join(latexFolderPath, "images/injAreaIViT.png"))

# plt.figure(figsize=(14, 9))
# plt.grid(True, which="both")
# plt.plot(t, var, lw=5, color='0.80', label="variance", marker="s", markersize=10)
# plt.plot(t, ivar, lw=5, color='0.40', label='inertial variance', marker="s", markersize=10)
# plt.plot(t, y, lw=5, color='blue', label='ref slope')
# plt.legend(loc="best", fontsize=30)