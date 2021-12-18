#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Eugenio Pescimoro
@email: eugenio.pescimoro@gmail.com
@Description: Inverse Gaussian CDF with two parameters
"""
###############################################################################
# The scipy.stats.invgauss.cdf is a standard scipy method to compute an Inverse Gaussian function
# however it only accepts one parameter (mu) because it assumes lambda = 1 -> we need 2 parameters 
#import scipy
#anaCdf = scipy.stats.invgauss.cdf(tt, mu1)
#axs1[0, 0].plot(tt, anaCdf)
#aa = scipy.stats.invgauss.fit(dCnorm)
###############################################################################
from scipy.stats import norm
from scipy import special
import numpy as np

def invGaussianPDF(t, l, mu, m):
    return m*np.sqrt(l/(2*np.pi*t*t*t))*np.exp(-(l*(t-mu)*(t-mu))/(2*mu*mu*t))

def invGaussianCDF(t, mu, l):
    return norm.cdf(np.sqrt(l/t)*(t/mu-1))+np.exp(2*l/mu)*norm.cdf(-np.sqrt(l/t)*(t/mu+1))

# From "van Genuchten M Th and Alves W J 1982 Analytical solu-tions of the one-dimensional convective-dispersive solute transport equation; US Dept. Agriculture Tech. Bull. No. 1661 151p." for X = X0
def invGauVG(u0, X, D0, T):
    return 0.5*special.erfc((X-u0*T)/(2*np.sqrt(D0*T)))+0.5*np.exp(u0*X/D0)*special.erfc((X+u0*T)/(2*np.sqrt(D0*T)))+0.5*(2+u0*X/D0+u0**2*T/D0)*np.exp(u0*X/D0)*special.erfc((X+u0*T)/(2*np.sqrt(D0*T)))-(u0**2*T/(np.pi*D0))**0.5*np.exp(u0*X/D0-(X+u0*T)**2/4*D0*T)