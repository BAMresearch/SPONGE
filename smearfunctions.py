# -*- coding: utf-8 -*-

import vtk
import pandas, scipy
import numpy as np
import scipy.spatial
import scipy.signal
import os, h5py

import scipy.stats # for the distributions
from scipy.interpolate import interp1d # for the interpolations
from scipy.stats import norm # gaussian distribution function
from .distfunctions import interpolate

def halfTrapzPDF(x, c, d):
    # this trapezoidal PDF is only defined from X >= 0, and is assumed
    # to be mirrored around that point.
    # Note that the integral of this PDF from X>0 will be 0.5.
    # source: van Dorp and Kotz, Metrika 2003, eq (1)
    # using a = -d, b = -c
    print("halfTrapzPDF called")
    assert(d > 0.)
    x = abs(x)
    pdf = x * 0.
    pdf[x < c] = 1.
    if d > c:
        pdf[(c <= x) & (x < d)] = (1./(d - c)) * (d - x[(c <= x) & (x < d)])
    norm = 1./(d + c)
    pdf *= norm
    return pdf, norm

def slitSmearTrapz(simDat, halfUmbra = 0.01, halfPenUmbra = 0.1, ndiv = 25):
    # this smears the profile according to slit smearing. It requires the input
    # of accurate beam profile parameters, AND it requires that the input simulated data
    # has a decent Guinier region at low-q (due to mirroring). 
    
    lims = [-halfPenUmbra, halfPenUmbra]
    x = np.linspace(lims[0], lims[1], ndiv)
    g, dummy = halfTrapzPDF(x, halfUmbra, halfPenUmbra)

    # for the numerical integration, we need to temporarily move away from the Pandas Dataframes
    I      = np.zeros([simDat.Q.values.size, g.size])
    IError = np.zeros([simDat.Q.values.size, g.size])

    for oi, offset in enumerate(x):
        dst = interpolate(simDat, QModifier = offset, modType = "slitSmear", mirror = True)
        
        # assign result into the to-be-integrated matrices
        I[:,oi] = dst.I.values * g[oi]
        IError[:,oi] = dst.IError.values * g[oi]
        
    print("Shape I: {}, shape x: {}".format(I.shape, x.shape))
    IInt = np.trapz(I, x = x[np.newaxis, :], axis = 1)
    IErrorInt = np.trapz(IError, x = x[np.newaxis, :], axis = 1)

    addDat = simDat.copy()
    addDat["Q"] = simDat["Q"] # reset Q
    addDat.I = IInt
    addDat.IError = IError
    
    return x, g, addDat
