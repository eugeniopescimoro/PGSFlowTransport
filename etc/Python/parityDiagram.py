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

FS = 1
sim = 1
interval = 1
U = [[] for i in range(sim)]
f = [[] for i in range(sim)] # Frequency or normalised probability
uf = [[] for i in range(sim)]
kf = [[] for i in range(sim)]
Kxx = [[] for i in range(sim)]
K = [[] for i in range(10)]
Ux = [[] for i in range(10)]
F = [[] for i in range(10)]

simPath = ['lowPeclet', 'mediumPeclet', 'highPeclet']

for i in range(0, sim, interval):
    # homeFolder = Path('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_3-highContrast/TS%d' % (FS+i))
    # homeFolder = Path('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_5-lowContrast/TS%d' % (FS+i))
    homeFolder = Path(os.path.join('/Users/pmxep5/Git/Hub/OpenFOAM/PGSFlowTransport/tutorials/Herten/varPeclet/', simPath[i]))
    latexFolderPath = Path('/Users/pmxep5/Git/Overleaf/Thesis')
    with open(Path(os.path.join(homeFolder, 'postProcessing/pdf/0/magUscaled-none_'))) as magUscaled:
        next(magUscaled)
        for index, line in enumerate(magUscaled):
            if float(line.split()[2])!=0:            
                U[i].append(float(line.split()[0]))
                Kxx[i].append(float(line.split()[1]))
                f[i].append(float(line.split()[2]))
    for index, fi in enumerate(f[i]):
        F[i].append(fi*U[i][index]*Kxx[i][index])

font = {'size': 30}
plt.rc('font', **font)
plt.figure(figsize=(14, 9)) # UNCOMMENT WHEN USING 3b OPTION
ax = plt.gca()
ax.scatter(U[0], Kxx[0], s=F[0])
plt.xlabel("$V^*_x$")
plt.ylabel("$K_x$")
ax.set_yscale('log')
ax.set_xscale('log')
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig(os.path.join(latexFolderPath, "images/parityDiagram.png"))

# for j in range(0, len(Ux)):
#     plt.loglog(Ux[j], F[j], color="%s" % col[j], label='Kxx=%s' % PermX[j])
# plt.axis([min(min(Ux, default=0), default=0), max(max(Ux, default=0), default=0), minYaxis, max(max(F, default=0), default=0)]) # np.compress = Pythonic way to slice list using boolean condition
# # plt.legend(loc="best")
# plt.xlabel("$V^*_x$")
# plt.ylabel("$p(V^*_x,K)$")
# plt.tight_layout()
# plt.grid(True, which="both")
# plt.savefig(os.path.join(latexFolderPath, "images/jojntPdfHerten.png"))