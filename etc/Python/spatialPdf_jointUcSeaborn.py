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
from itertools import islice

FS = 1
sim = 1
interval = 1
nClass = 50 # Number of classes for the fields (k, U, c) to be divided into (e.g. lowContrast=4, highContrast=4, joint(K, U) for Herten=10, joint(c, U) for Herten=4)
c = [[] for i in range(sim)]
U = [[] for i in range(sim)]
f = [[] for i in range(sim)] # Frequency or normalised probability
uf = [[] for i in range(sim)]
kf = [[] for i in range(sim)]
Kxx = [[] for i in range(sim)]
K = [[] for i in range(nClass)]
Ux = [[] for i in range(nClass)]
F = [[] for i in range(nClass)]
C = [[] for i in range(nClass)]

font = {'size': 30}
plt.rc('font', **font)
# plt.figure(figsize=(14, 9)) # UNCOMMENT WHEN USING 3b OPTION
simPath = ['lowPeclet', 'mediumPeclet', 'highPeclet']

subprocess.run(['/bin/bash', '-c', '. /opt/openfoam8/etc/bashrc']) # It loads the environment variable OpenFOAM
for i in range(0, sim, interval):
# 0) Set the path for the bash script to be executed, usually the testcase home folder
    # homeFolder = Path('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_3-highContrast/TS%d' % (FS+i))
    # homeFolder = Path('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_5-lowContrast/TS%d' % (FS+i))
    homeFolder = Path(os.path.join('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/Herten/varPeclet/', simPath[i]))
    latexFolderPath = Path('/home/pmxep5/Git/Overleaf/Thesis')

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
#     Kxx
#     {
#         // Load the library containing the 'coded' functionObject
#         functionObjectLibs ("libutilityFunctionObjects.so");
        
#         type coded;        
#         // Name of on-the-fly generated functionObject
#         name symmTensorXXcomponent;
        
#         redirectType Kxx;
        
#         codeExecute
#         #{
#             Info<< "Looking up filed K\n" << endl;
#             // Lookup K
#             const volSymmTensorField& K = mesh().lookupObject<volSymmTensorField>("K");
#             volScalarField Kx("Kx", K.component(tensor::XX));
#             Kx.write();
#         #};
#     }
#     c
#     {
#         // Load the library containing the 'coded' functionObject
#         functionObjectLibs ("libutilityFunctionObjects.so");
        
#         type coded;        
#         // Name of on-the-fly generated functionObject
#         name symmTensorXXcomponent;
        
#         redirectType c;
        
#         codeExecute
#         #{
#             Info<< "Looking up filed c\n" << endl;
#             // Lookup c
#             const volScalarField& c = mesh().lookupObject<volScalarField>("c");
#             //volScalarField Kx("Kx", K.component(tensor::XX));
#             //c.write();
#         #};
#     }
# }

# 2a) Launch postProcess
    # os.chdir(homeFolder)
    # functionObject and spatialPdf yield the PDF of the velocity spatial field
    # subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 postProcess -time 0 -field U -parallel -func magScaled']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
    # subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 spatialPdf -time 0 -field magUscaled -parallel -logbin -nbins 50']) # It runs the spatialPdf post processing utility and stores the output in postProcessing/png/0/magUscaled-none_
    # functionObject and spatialPdf are run to obtain the conditional PDF of the spatial velocity fields given a facies permeability
    # subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 postProcess -time 0 -fields \'(U K)\' -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/Kx
    # subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 spatialPdf -time 0 -field magUscaled -parallel -logbin -nbins 50 -joint Kx']) # It runs the spatialPdf post processing utility with a weight function and stores the output in postProcessing/pdf/0/magUscaled-Kx_
    
# # 2b) Plot the output with gnuplot: make sure the plotPdf gnuplot filed exist otherwise create one and copy and paste the following instructions
# set logscale x
# set logscale y
# set title "U/Uave distribution"
# set xlabel 'U/Uave [m/s]'
# set ylabel 'Normalised probability'
# plot "< tail -n +2 postProcessing/pdf/0/magUscaled-none_" using 1:3 title "H Vel Dist" with lines
# pause 1
# reread
# # subprocess.run(['/bin/bash', '-c', 'gnuplot plotPdf'])

# # 3a) Plot spatial pdf with Python
#     os.chdir(homeFolder)
#     # N.M.FUCKING B.: FOR SOME UNKNOWN REASON IF NUMBER OF "PHYSICAL PROCESSOR < mpirun -np" THE postProcessing AND spatialPdf FUNCTIONS NEED TO BE RUN MANUALLY FROM THE TERMINAL!!
#     # functionObject and spatialPdf yield the PDF of the velocity spatial field
#     subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 postProcess -func magScaled -dict system/fieldMetricsDict -field U -time 0 -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
#     subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 0 -logbin -nbins 50']) # It runs the spatialPdf post processing utility and stores the output in postProcessing/png/0/magUscaled-none_
#     col = ['0.0', '0.2', '0.4', '0.6']
#     lab = ['0.4', '0.6', '0.8', '1.0']
#     minYaxis = 1e-3
#     with open(Path(os.path.join(homeFolder, 'postProcessing/pdf/0/magUscaled-none_'))) as magUscaled:
#         next(magUscaled)
#         for index, line in enumerate(magUscaled):
#             U[i].append(float(line.split()[0]))
#             f[i].append(float(line.split()[2]))
#             uf[i].append(U[i][index]*f[i][index])
#     plt.loglog(U[i], uf[i], color="%s" % col[i], label='Lx=%s' % lab[i])
#     plt.axis([np.compress(np.array(uf[i])>minYaxis, U[i])[0], max(U[i]), minYaxis, max(uf[i])]) # np.compress = Pythonic way to slice list using boolean condition
#     # plt.legend(loc="best")
#     plt.xlabel("$V^*_x$")
#     plt.ylabel("$p(V^*_x)$")
#     plt.grid(True, which="both")
#     plt.tight_layout()
#     # plt.savefig(os.path.join(latexFolderPath, "images/magUscaledHerten.png"))
    
# 3b) Plot joint spatial pdf with Python
    os.chdir(homeFolder)
    # functionObject and spatialPdf are run to obtain the conditional PDF of the spatial velocity fields given a facies permeability or the concentration field
    # FOR -np > number  of physical processors, the following command needs to be run manually: mpirun --oversubscribe -np 96 postProcess -funcs '(magScaled c)' -dict system/fieldMetricsDict -fields '(U c)' -time 100000 -parallel
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 postProcess -funcs \'(magScaled c)\' -dict system/fieldMetricsDict -fields \'(U c)\' -time 100000 -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 postProcess -funcs \'(magScaled Kxx)\' -dict system/fieldMetricsDict -fields \'(U K)\' -time 0 -parallel']) # It runs the functionObject added in 1) and the output is stored as Kx and magUscaled in processor*/0/
    # FOR -np > number  of physical processors, the following command needs to be run manually: mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 100000 -logbin -nbins 50 -joint c
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 100000 -logbin -nbins 50 -joint c'])
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 0 -logbin -nbins 50 -joint Kx']) # It runs the spatialPdf post processing utility with a weight function and stores the output in postProcessing/pdf/0/magUscaled-Kx_

    with open(Path(os.path.join(homeFolder, 'postProcessing/pdf/100000/magUscaled-none_'))) as magUscaled:
        next(magUscaled)
        for index, line in enumerate(magUscaled):
            # if float(line.split()[2])!=0:            
                U[i].append(float(line.split()[0]))
                c[i].append(float(line.split()[1]))
                f[i].append(float(line.split()[2]))

#         for index, line in enumerate(magUscaled):
#             f[i].append(float(line.split()[2]))                
# length_to_split = [50 for i in range(nClass)]
# Inputt = iter(f[i])
# Output = [list(islice(Inputt, elem)) for elem in length_to_split]

# # Bins for joint distribution (V, c)
#     # cClasses = np.linspace(min(c[0]), max(c[0]), num=nClass+1)
#     cClasses = np.logspace(math.log10(min(c[0])), math.log10(max(c[0])), num=nClass+1)
#     for index, x in enumerate(c[i]):
#         for j in range(nClass):
#             if cClasses[j] < x <= cClasses[j+1]:
#                 Ux[j].append(U[i][index])
#                 C[j].append(x)
#                 F[j].append(f[i][index]*Ux[j][-1]*x)                
# Uave = [np.mean(Ui) for Ui in Ux]
# Cave = [np.mean(ci) for ci in C]
# Fi = [sum(Fi) for Fi in F]

import seaborn as sns
# sns.heatmap(Output)
# UcJointPdf = sns.jointplot(Uave, Cave, Fi, kind="hist")
UcJointPdf = sns.jointplot(U[0], c[0], f[0], kind="hist")
UcJointPdf.ax_joint.set_xscale('log')
UcJointPdf.ax_joint.set_yscale('log')
plt.grid(True, which="both")
UcJointPdf.set_axis_labels('U*', 'c', fontsize=30)
plt.savefig(os.path.join(latexFolderPath, "images/jointUcColourMap.png"))