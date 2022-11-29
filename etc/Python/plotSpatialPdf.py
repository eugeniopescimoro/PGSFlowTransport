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

font = {'size': 50}
plt.rc('font', **font)
plt.figure(figsize=(14, 9)) # UNCOMMENT WHEN USING 3b OPTION
simPath = ['lowPeclet', 'mediumPeclet', 'highPeclet']

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

# 3a) Plot marginal spatial pdf with Python
    # os.chdir(homeFolder)
    # # functionObject and spatialPdf yield the PDF of the velocity spatial field
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 postProcess -dict system/fieldMetricsDict -field U -time 0 -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
    # # subprocess.run(['/bin/bash', '-c', 'mpirun -np 96 postProcess -time 0 -field U -parallel -func magScaled']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
    # # mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 0 -logbin -nbins 50
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 0 -logbin -nbins 50'])
    # # subprocess.run(['/bin/bash', '-c', 'mpirun -np 96 spatialPdf -time 0 -field magUscaled -parallel -logbin -nbins 50']) # It runs the spatialPdf post processing utility and stores the output in postProcessing/png/0/magUscaled-none_
    # # mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 0 -logbin -nbins 50
    # col = ['0.0', '0.2', '0.4', '0.6']
    # lab = ['0.4', '0.6', '0.8', '1.0']
    # minYaxis = 1e-3
    # with open(Path(os.path.join(homeFolder, 'postProcessing/pdf/0/magUscaled-none_'))) as magUscaled:
    #     next(magUscaled)
    #     for index, line in enumerate(magUscaled):
    #         U[i].append(float(line.split()[0]))
    #         f[i].append(float(line.split()[2]))
    #         uf[i].append(U[i][index]*f[i][index])
    # plt.loglog(U[i], uf[i], color="%s" % col[i], label='Lx=%s' % lab[i])
    # plt.axis([np.compress(np.array(uf[i])>minYaxis, U[i])[0], max(U[i]), minYaxis, max(uf[i])]) # np.compress = Pythonic way to slice list using boolean condition
    # # plt.legend(loc="best")
    # plt.xlabel("$V^*_x$")
    # plt.ylabel("$p(V^*_x)$")
    # plt.grid(True, which="both")
    # plt.tight_layout()
    # plt.savefig(os.path.join(latexFolderPath, "images/magUscaledHerten.png"))
    
# # 3b) Plot joint spatial pdf with Python
    os.chdir(homeFolder)
    # functionObject and spatialPdf are run to obtain the conditional PDF of the spatial velocity fields given a facies permeability
    subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 postProcess -dict system/fieldMetricsDict -fields \'(U c)\' -time 100000 -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 postProcess -dict system/fieldMetricsDict -fields \'(U K)\' -time 0 -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/magUscaled
    # subprocess.run(['/bin/bash', '-c', 'mpirun -np 96 postProcess -time 0 -fields \'(U K)\' -parallel']) # It runs the functionObject added in 1) and the output is stored in processor*/0/Kx
    subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 100000 -logbin -nbins 50 -joint c'])
    # subprocess.run(['/bin/bash', '-c', 'mpirun --oversubscribe -np 96 spatialPdf -parallel -field magUscaled -time 0 -logbin -nbins 50 -joint Kx'])
    # subprocess.run(['/bin/bash', '-c', 'mpirun -np 96 spatialPdf -time 0 -field magUscaled -parallel -logbin -nbins 50 -joint Kx']) # It runs the spatialPdf post processing utility with a weight function and stores the output in postProcessing/pdf/0/magUscaled-Kx_
    # NB: THE FIRST 4 "IF" ARE THOUGHT FOR HIGH PERMEABILITY CONTRAST WHILE THE SECOND 4 "IF" ARE NEEDED FOR LOW CONTRAST PERMEABILITY 
    # Labels for the paper
    # PermX = ['1e-10', '1e-11', '1e-12', '1e-13']    
    # PermX = ['1e-9', '1e-11', '1e-13', '1e-15']
    # Labels for Herten
    PermX = ['1.18e-8', '8.62e-9', '2.36e-9', '2.09e-9', '2.09e-10', '2.27e-11', '2.09e-11', '1.27e-11', '5.53e-12', '3.90e-12', '5.44e-14']
    # col = ['0.0', '0.2', '0.4', '0.6']
    col = ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9']
    minYaxis = 1e-3
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
                                        
############################################################   
# Bins for high contrast         
        # if x >= 1e-9:
        #     K[0].append(x)
        #     Ux[0].append(U[i][index])
        #     F[0].append(f[i][index]*Ux[0][-1]*x)
        # if 1e-10 >= x >= 1e-12:
        #     K[1].append(x)
        #     Ux[1].append(U[i][index])
        #     F[1].append(f[i][index]*Ux[1][-1]*x)
        # if 1e-12 >= x >= 1e-14:
        #     K[2].append(x)
        #     Ux[2].append(U[i][index])
        #     F[2].append(f[i][index]*Ux[2][-1]*x)
        # if x <= 1e-14:
        #     K[3].append(x)
        #     Ux[3].append(U[i][index])
        #     F[3].append(f[i][index]*Ux[3][-1]*x)
############################################################
# Bins for low contrast
        # if x >= 1e-10:
        #     K[0].append(x)
        #     Ux[0].append(U[i][index])
        #     F[0].append(f[i][index]*Ux[0][-1]*x)
        # if 1e-10 > x >= 1e-11:
        #     K[1].append(x)
        #     Ux[1].append(U[i][index])
        #     F[1].append(f[i][index]*Ux[1][-1]*x)
        # if 1e-11 > x >= 1e-12:
        #     K[2].append(x)
        #     Ux[2].append(U[i][index])
        #     F[2].append(f[i][index]*Ux[2][-1]*x)
        # if x < 1e-12:
        #     K[3].append(x)
        #     Ux[3].append(U[i][index])
        #     F[3].append(f[i][index]*Ux[3][-1]*x)        
for j in range(0, len(Ux)):
    plt.loglog(Ux[j], F[j], color="%s" % col[j], label='Kxx=%s' % PermX[j])
plt.axis([min(min(Ux, default=0), default=0), max(max(Ux, default=0), default=0), minYaxis, max(max(F, default=0), default=0)]) # np.compress = Pythonic way to slice list using boolean condition
# plt.legend(loc="best")
plt.xlabel("$V^*_x$")
plt.ylabel("$p(V^*_x,K)$")
plt.tight_layout()
plt.grid(True, which="both")
# plt.savefig(os.path.join(latexFolderPath, "images/jointPdfHerten.png"))

###############################################################################
# import seaborn as sns
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
    # UKjointPdf.savefig(os.path.join(latexFolderPath, "images/jointUKpdfHC.png"))