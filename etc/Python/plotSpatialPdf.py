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

FS = 1
sim = 4
interval = 1
U = [[] for i in range(sim)]
f = [[] for i in range(sim)]
# U = []
# f = []

font = {'size': 30}
plt.rc('font', **font)
plt.figure(figsize=(14, 9))

for i in range(0, sim, interval):
# 0) Set the path for the bash script to be executed, usually the testcase home folder
    homeFolder = Path('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_3-highContrast/TS%d' % (FS+i))
    latexFolderPath = Path('/home/pmxep5/OneDrive/Nottingham/Write/Articles/PGSFoam/')

# # 1) Verify the following functionObject to be present in system/controlDict
# functions
# {
#     magScaled
#     {
#         // Load the library containing the 'coded' functionObject
#         functionObjectLibs ("libutilityFunctionObjects.so");
        
#         type coded;        
#         // Name of on-the-fly generated functionObject
#         name UmagScaled;
        
#         redirectType magScaled;
        
#         codeExecute
#         #{
#             Info<< "Looking up filed U\n" << endl;
#             // Lookup U
#             const volVectorField& U = mesh().lookupObject<volVectorField>("U"); 
#             scalar vol(gSum(mesh().V()));
#             volScalarField magU(mag(U));
#             volScalarField magUscaled("magUscaled", magU/gSum(magU.primitiveField()*mesh().V())*vol);
#             magUscaled.write();
#         #};
#     }
# }

# 2) Launch postProcess
    os.chdir(homeFolder)
    subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 postProcess -time 0 -field U -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
    subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 spatialPdf -time 0 -field magUscaled -parallel -logbin -nbins 50']) # It runs the spatialPdf post processing utility and stores the output in postProcessing/pdf/0/magUscaled-none_

# # 3a) Plot the output with gnuplot: make sure the plotPdf gnuplot filed exist otherwise create one and copy and paste the following instructions
# set logscale x
# set logscale y
# set title "U/Uave distribution"
# set xlabel 'U/Uave [m/s]'
# set ylabel 'Normalised probability'
# plot "< tail -n +2 postProcessing/pdf/0/magUscaled-none_" using 1:3 title "H Vel Dist" with lines
# pause 1
# reread
# # subprocess.run(['/bin/bash', '-c', 'gnuplot plotPdf'])

# 3b) Plot spatial pdf with Python
    length = ['0.4', '0.6', '0.8', '1.0']
    with open(Path(os.path.join(homeFolder, 'postProcessing/pdf/0/magUscaled-none_'))) as magUscaled:
        next(magUscaled)
        for line in magUscaled:
            U[i].append(float(line.split()[0]))
            f[i].append(float(line.split()[2]))
    plt.loglog(U[i], f[i], label='Lx=%s' % length[i])
    plt.legend(loc="best")
    plt.xlabel('Vx/Vave')
    plt.ylabel('Normalised probability')
    # plt.savefig(os.path.join(latexFolderPath, "images/magUscaled.pdf"))