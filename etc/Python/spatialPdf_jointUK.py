#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Parse and plot the pdf of the velocity field computed by spatialPdf
"""
from pathlib import Path
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import math

FS = 1
sim = 1
interval = 1
nClass = 10 # Number of classes for the fields (k, U, c) to be divided into (e.g. lowContrast=4, highContrast=4, Herten=10)
U = [[] for i in range(sim)]
f = [[] for i in range(sim)] # Frequency or normalised probability
uf = [[] for i in range(sim)]
kf = [[] for i in range(sim)]
Kxx = [[] for i in range(sim)]
K = [[] for i in range(nClass)]
Ux = [[] for i in range(nClass)]
F = [[] for i in range(nClass)]
c = [[] for i in range(nClass)]

font = {'size': 30}
plt.rc('font', **font)
plt.figure(figsize=(14, 9)) # UNCOMMENT WHEN USING 3b OPTION
simPath = ['lowPeclet', 'mediumPeclet', 'highPeclet']

subprocess.run(['/bin/bash', '-c', '. /opt/openfoam8/etc/bashrc']) # It loads the environment variable OpenFOAM
for i in range(0, sim, interval):
# 0) Set the path for the bash script to be executed, usually the testcase home folder
    # homeFolder = Path('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_3-highContrast/TS%d' % (FS+i))
    # homeFolder = Path('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_5-lowContrast/TS%d' % (FS+i))
    homeFolder = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/Herten/varPeclet/', simPath[i]))
    latexFolderPath = Path('/home/pmxep5/Git/Overleaf/Thesis')

# N.M.FUCKING B.: FOR SOME UNKNOWN REASON IF NUMBER OF "PHYSICAL PROCESSOR < mpirun -np" THE postProcessing AND spatialPdf FUNCTIONS NEED TO BE RUN MANUALLY FROM THE TERMINAL!!
    
# 3b) Plot joint spatial pdf with Python
    os.chdir(homeFolder)
    # functionObject and spatialPdf are run to obtain the conditional PDF of the spatial velocity fields given a facies permeability
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 postProcess -funcs \'(magScaled c)\' -dict system/fieldMetricsDict -fields \'(U c)\' -time 100000 -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
    subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 postProcess -funcs \'(magScaled Kxx)\' -dict system/fieldMetricsDict -fields \'(U K)\' -time 0 -parallel']) # It runs the functionObject added in 1) and the output is stored as Kx and magUscaled in processor*/0/
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 100000 -logbin -nbins 50 -joint c'])
    subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 0 -logbin -nbins 50 -joint Kx']) # It runs the spatialPdf post processing utility with a weight function and stores the output in postProcessing/pdf/0/magUscaled-Kx_
    # Labels for Herten
    PermX = ['1.18e-8', '8.62e-9', '2.36e-9', '2.09e-9', '2.09e-10', '2.27e-11', '2.09e-11', '1.27e-11', '5.53e-12', '3.90e-12', '5.44e-14']
    # col = ['0.0', '0.2', '0.4', '0.6']
    col = ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']
    minYaxis = 1e-2
    with open(Path(os.path.join(homeFolder, 'postProcessing/pdf/0/magUscaled-none_'))) as magUscaled:
        next(magUscaled)
        for index, line in enumerate(magUscaled):
            if float(line.split()[2])!=0:            
                U[i].append(float(line.split()[0]))
                Kxx[i].append(float(line.split()[1]))
                f[i].append(float(line.split()[2]))
###########################################################
# Bins for Herten     
    for index, x in enumerate(Kxx[i]):
        if 1.18e-8 >= x > 8.62e-9:
            K[0].append(x)
            Ux[0].append(U[i][index])
            F[0].append(f[i][index]*Ux[0][-1]*x)
        else:
            if 8.62e-9 >= x > 2.36e-9:
                K[1].append(x)
                Ux[1].append(U[i][index])
                F[1].append(f[i][index]*Ux[1][-1]*x)
            else:
                if 2.36e-9 >= x > 2.09e-10:
                    K[2].append(x)
                    Ux[2].append(U[i][index])
                    F[2].append(f[i][index]*Ux[2][-1]*x)
                else:
                    if 2.09e-10 >= x > 2.27e-11:
                        K[3].append(x)
                        Ux[3].append(U[i][index])
                        F[3].append(f[i][index]*Ux[3][-1]*x)
                    else:
                        if 2.27e-11 >= x > 2.09e-11:
                            K[4].append(x)
                            Ux[4].append(U[i][index])
                            F[4].append(f[i][index]*Ux[4][-1]*x)
                        else:
                            if 2.09e-11 >= x > 1.27e-11:
                                K[5].append(x)
                                Ux[5].append(U[i][index])
                                F[5].append(f[i][index]*Ux[5][-1]*x)
                            else:
                                if 1.27e-11 >= x > 5.53e-12:
                                    K[6].append(x)
                                    Ux[6].append(U[i][index])
                                    F[6].append(f[i][index]*Ux[6][-1]*x)
                                else:
                                    if 5.53e-12 >= x > 3.90e-12:
                                        K[7].append(x)
                                        Ux[7].append(U[i][index])
                                        F[7].append(f[i][index]*Ux[7][-1]*x)
                                    else:
                                        if 3.90e-12 >= x > 5.44e-14:
                                            K[8].append(x)
                                            Ux[8].append(U[i][index])
                                            F[8].append(f[i][index]*Ux[8][-1]*x)
                                        else:
                                            if 5.44e-14 >= x:
                                                K[9].append(x)
                                                Ux[9].append(U[i][index])
                                                F[9].append(f[i][index]*Ux[9][-1]*x)                                        
###########################################################                               
# Plot joint (k, U)
for j in range(0, len(Ux)):
    Ux[j] = Ux[j] or [0] # If list is empty it sets it to 0 otherwise error in the min/max function
    F[j] = F[j] or [0] # If list is empty it sets it to 0 otherwise error in the min/max function
    plt.loglog(Ux[j], F[j], color="%s" % col[j], label='Kxx=%s' % PermX[j])
plt.axis([min(min(Ux, key=min)), max(max(Ux, key=max)), minYaxis, max(max(F, key=max))]) # np.compress = Pythonic way to slice list using boolean condition
# plt.legend(loc="best")
plt.xlabel("$V^*_x$")
plt.ylabel("$p(V^*_x, Kxx)$")
plt.tight_layout()
plt.grid(True, which="both")
# plt.savefig(os.path.join(latexFolderPath, "images/jointPdfHerten.png"))