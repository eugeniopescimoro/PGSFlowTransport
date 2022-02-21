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
import seaborn as sns
import numpy as np

FS = 2
sim = 1
interval = 1
U = [[] for i in range(sim)]
f = [[] for i in range(sim)] # Frequency or normalised probability
uf = [[] for i in range(sim)]
kf = [[] for i in range(sim)]
Kxx = [[] for i in range(sim)]
K = [[] for i in range(4)]
Ux = [[] for i in range(4)]
F = [[] for i in range(4)]

font = {'size': 30}
plt.rc('font', **font)
plt.figure(figsize=(14, 9)) # UNCOMMENT WHEN USING 3b OPTION

for i in range(0, sim, interval):
# 0) Set the path for the bash script to be executed, usually the testcase home folder
    homeFolder = Path('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_3-highContrast/TS%d' % (FS+i))
    # homeFolder = Path('/data/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp/scat_5-lowContrast/TS%d' % (FS+i))
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
# }

# 2) Launch postProcess
    os.chdir(homeFolder)
    # subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 postProcess -time 0 -field U -parallel -func magScaled']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
    subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 postProcess -time 0 -fields \'(U K)\' -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/Kx
    # subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 spatialPdf -time 0 -field magUscaled -parallel -logbin -nbins 50']) # It runs the spatialPdf post processing utility and stores the output in postProcessing/pdf/0/magUscaled-none_
    subprocess.run(['/bin/bash', '-c', 'mpirun -np 12 spatialPdf -time 0 -field magUscaled -parallel -logbin -nbins 50 -joint Kx']) # It runs the spatialPdf post processing utility with a weight function and stores the output in postProcessing/pdf/0/magUscaled-Kx_
    
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
    # length = ['0.4', '0.6', '0.8', '1.0']
    # with open(Path(os.path.join(homeFolder, 'postProcessing/pdf/0/magUscaled-none_'))) as magUscaled:
    #     next(magUscaled)
    #     for index, line in enumerate(magUscaled):
    #         U[i].append(float(line.split()[0]))
    #         f[i].append(float(line.split()[2]))
    #         uf[i].append(U[i][index]*f[i][index])
    # plt.loglog(U[i], uf[i], label='Lx=%s' % length[i])
    # plt.axis([min(U[i]), max(U[i]), min(uf[i]), max(uf[i])])
    # plt.legend(loc="best")
    # plt.xlabel('Vx/Vave')
    # plt.ylabel('Normalised joint probability*Vx/Vave')
    # # plt.savefig(os.path.join(latexFolderPath, "images/magUscaledLC.pdf"))
    
# 3c) Plot joint spatial pdf with Python
    PermX = ['1e-9', '1e-11', '1e-13', '1e-15']
    with open(Path(os.path.join(homeFolder, 'postProcessing/pdf/0/magUscaled-none_'))) as magUscaled:
        next(magUscaled)
        for index, line in enumerate(magUscaled):
            if float(line.split()[2])!=0:            
                U[i].append(float(line.split()[0]))
                Kxx[i].append(float(line.split()[1]))
                f[i].append(float(line.split()[2]))
    
    for index, x in enumerate(Kxx[i]):
        if x >= 1e-9:
            K[0].append(x)
            Ux[0].append(U[i][index])
            F[0].append(f[i][index]*Ux[0][-1]*x)
        if 1e-10 >= x >= 1e-12:
            K[1].append(x)
            Ux[1].append(U[i][index])
            F[1].append(f[i][index]*Ux[1][-1]*x)
        if 1e-12 >= x >= 1e-14:
            K[2].append(x)
            Ux[2].append(U[i][index])
            F[2].append(f[i][index]*Ux[2][-1]*x)
        if x <= 1e-14:
            K[3].append(x)
            Ux[3].append(U[i][index])
            F[3].append(f[i][index]*Ux[3][-1]*x)
    for j in range(0, len(Ux)):
        plt.loglog(Ux[j], F[j], label='Kxx=%s' % PermX[j])
    plt.legend(loc="best")
    plt.xlabel('Vx/Vave')
    plt.ylabel('Norm joint prob*Vx/Vave*Kxx')
    plt.savefig(os.path.join(latexFolderPath, "images/jointPdfdHC.pdf"))

    # for j in range (0, len(U[i])):       
    #     uf[i].append(f[i][j]*U[i][j]*Kxx[i][j])
    # UKjointPdf = sns.jointplot(U[i], Kxx[i], uf[i], kind="hist")

    # for j in range (0, len(U[i])):       
    #     uf[i].append(U[i][j]*f[i][j])
    #     kf[i].append(Kxx[i][j]*f[i][j])
    # UKjointPdf = sns.jointplot(uf[i], kf[i], f[i], kind="hist")
    
    # UKjointPdf = sns.jointplot(U[i], Kxx[i], f[i], kind="hist")
    # UKjointPdf = sns.jointplot(np.log10(U[i]), np.log10(Kxx[i]), np.log10(f[i]), kind="hist")
    # UKjointPdf.set_axis_labels('log(U/Uave)', 'log(Kxx)', fontsize=24)
    # sns.jointplot(U[i], Kxx[i], f[i], kind="hist")
    # UKjointPdf.savefig(os.path.join(latexFolderPath, "images/jointUKpdfHC.pdf"))







