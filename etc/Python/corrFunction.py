#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Mar 12 11:09:26 2021

@author: Eugenio Pescimoro
email: eugenio.pescimoro@gmail.com
Plot (auto)correlation function as output from OpenFOAM fieldMetrics functionObject 
"""

# Import section ############################################################## 
import os
from matplotlib import pyplot as plt
#import numpy as np 
from matplotlib.ticker import MultipleLocator

# Draw blank canvas ###########################################################
font = {'size': 24}
plt.rc('font', **font)
corrFun = plt.figure(figsize=(14, 9))
axes = corrFun.add_axes([0.15, 0.1, 0.8, 0.8]) # [left, bottom, width, height]
axes.set_xlabel("Dist [m]")
axes.set_ylabel("Correlation value [m]")
plt.title(label = "(Auto)correlation function")
axes.xaxis.set_major_locator(MultipleLocator(0.1))
axes.tick_params(which='major', length=7, width=2)
axes.xaxis.set_major_formatter('{x:1.1f}')
axes.xaxis.set_minor_locator(MultipleLocator(0.01))
axes.xaxis.grid(True, which='both')

# Loop through simulation folders and parse autocorrelation_K file ############
for i in range(1, 2, 1):
    os.chdir('/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_3/TS%d/postProcessing/fieldMetrics/0' % i)
    x = []
    corr = []
    with open('autocorrelation_K') as corrK:
        for line in corrK: # It parses all line in autocorrelation_K
            x.append(float(line.split()[0]))
            corr.append(float(line.split()[1]))
    plt.plot(x, corr, label='TS%d' % i)

# Add legend and save image ###################################################
plt.legend(loc = 'best')
corrFun.savefig("/home/pmxep5/OpenFOAM/pmxep5-8/PGSFlowTransport/tutorials/RESULTS/stopConcAdapTmstp_3/images/correlationFuncK.pdf")

#plt.xlabel("Dist [m]")
#plt.ylabel("Correlation value [-]")
#plt.title('(Auto)correlation function')
#plt.show()
#plt.savefig("/home/pmxep5/OneDrive/Nottingham/Results/Images/veLcorrTransport/correlationFuncK.pdf")
