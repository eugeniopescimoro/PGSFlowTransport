#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Parse script for transport properties and log OpenFOAM file (setRandomField, simpleDarcyFoam, scalarTransportFoam)
"""
import re
import os
import numpy as np

def parseLog(logPath, s, cl, dd, mvel, c, t, m, kvol, kval):
    mass = []
    time = []
    conc = []
    os.chdir(logPath)          
    with open("log") as log:
        parse = False
        oneVel = True
        # for idx, line in enumerate(log): # It reads all the lines in the log file and it filters out only those for which idx % 100 is not 0 (reminder of idx/100 is not 0)        
        for line in log:
            #if "scalarTransportFoam" in line:
            if any(s in line for s in ("adaptiveScalarTransportFoam", "macroScalarTransport")):
                parse = True
            if "postProcess" in line:
                parse = False
            if "volumes:" in line:
                kvol.append([float(i) for i in re.split(' |\(|\)', line)[slice(2, 6)]])
            if "field values:" in line:
                kval.append([float(i) for i in re.split(' |\(|\)', line)[slice(3, 7)]])
            if parse == True and "Total mass =" in line:
                mass.append(float(line.split()[3]))
            if parse == True and oneVel == True and "Mean vel =" in line:
                mvel.append([float(i) for i in re.split(' |\(|\)', line)[slice(4, 7)]]) # It does not "append" the values since steady state flow has constant velocities for all time steps
                # NMB: re.split COUNTS THE WHITE SPACES WHILE line.split NO!
                oneVel = False
            if parse == True and "Flux out =" in line:
                conc.append(float(line.split()[3]))
            if parse == True and "Time =" in line:
                time.append(float(line.split()[2]))
            if parse == True and "Adaptive time =" in line:
                time.append(float(line.split()[3]))
            if "boundingBox:" in line:
                dd.append([float(i) for i in re.split(' |\(|\)', line)[slice(9, 12)]])
            if re.search(r'\b%s\b' % (re.escape("Lcorr")), line): # Search for the WHOLE word "Lcorr" and not for the occurrencies within other string, e.g. veLcorrTransport                
                cl.append([float(i) for i in re.split(' |\(|\)', line)[slice(2, 5)]])
        # Select one element of the list every 100 to reduce the memory usage  

        # Select concentrations above threshold and correspondent times
        # cBoolean = np.logical_and(np.array(conc)>1e-6, np.array(conc)<1) 
        # c = [val for i, val in enumerate(conc) if cBoolean[i]]
        # t = [val for i, val in enumerate(time) if cBoolean[i]] # it selects the time only if cBoolean is True
        
        if len(time)>len(conc):
            time = time[:-1]        
        mass = mass[0:-1:s]
        m.append(np.array(mass))
        time = time[0:-1:s]
        t.append(np.array(time))
        conc = conc[0:-1:s]
        c.append(np.array(conc))

def parseConstants(transPropPath):
    D = [] 
    mu = []
    rho = []
    g = []
    os.chdir(transPropPath)
    with open("constant/transportProperties") as TP:
        for line in TP:
            if any(s in line for s in ("DT", "Dm")):
                D.append(float(str(line.split()[9]).replace(';','')))
            if "mu" in line:
                mu.append(float(str(line.split()[9]).replace(';','')))
            if "rho" in line:
                rho.append(float(str(line.split()[9]).replace(';','')))
    with open("constant/g") as G:
        for line in G:
            if "value" in line:
                g = [float(i) for i in re.split(' |\(|\)', line)[slice(12, 15)]]
    return D, mu, rho, g

def parseInitialConditions(inCondPath):
    Pin = [] 
    Pout = []
    os.chdir(inCondPath)
    with open("0/p.orig") as Porig:
        inlet = False
        outlet = False
        for line in Porig:
            if "inlet" in line:
                inlet = True
                outlet = False
            if "outlet" in line:
                inlet = False
                outlet = True
            if inlet == True and "value uniform" in line:
                Pin.append(float(str(line.split()[2]).replace(';','')))
            if outlet == True and "value uniform" in line:
                Pout.append(float(str(line.split()[2]).replace(';','')))
    deltaPx = Pout[0]-Pin[0]
    return deltaPx

def parseSetFieldsDict(setFieldPath):
    Xbox = [] # Variable which stores the horizontal distance at which the concentration source is placed
    os.chdir(setFieldPath)
    parse = False
    with open("system/topoSetDict") as SFD:
        for line in SFD:
            if "boxToCell" in line: # To discard the occurrencies of the word "box" in the comments
                parse = True            
            if parse == True and re.search(r'\b%s\b' % (re.escape("box")), line): # Search for the WHOLE word "box" and not for the occurrencies within other string, e.g. boxToCell
                Xbox1 = [float(i) for i in re.split(' |\(|\)', line)[slice(14, 15)]]
                Xbox2 = [float(i) for i in re.split(' |\(|\)', line)[slice(19, 20)]]
    Xbox = (float(Xbox1[0])+float(Xbox2[0]))/2
    return Xbox
    