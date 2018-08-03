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

def distSpreadGaussian(simDat, sigma = 0.01, ndiv = 10):
    # we program a dilation-/contraction-type polydispersity in here.
    lims = scipy.stats.distributions.norm.ppf([0.05, 0.95], loc = 1, scale = sigma)
    x = np.linspace(lims[0], lims[1], ndiv)
    g = scipy.stats.distributions.norm.pdf(x, loc = 1, scale = sigma)
    addDat = simDat * 0 # initialize as zero
    for oi, offset in enumerate(x):
        dst = interpolate(simDat, QModifier = offset, modType = "multiply")
        addDat += dst * g[oi] # / g.sum()
        addDat["Q"] = simDat["Q"] 
    return x, g, addDat

def interpolate(dataset, QModifier = 1., modType = "multiply", mirror = False): 
    # function modified from imp2
    # helper function to take care of interpolating the intensity within the limits
    # can do scaling (dilation, contraction) for size distributions, or shifting for
    # beam smearing calculations. For the latter, it is recommended to mirror the 
    # intensity to negative Q values, using mirror = True. This only works for 
    # intensity profiles that are reasonably flat at q < qmin!
    # A third option is modType = "slitSmear", which calculates using 
    # Q = np.sqrt(Q**2 + qModifier**2)
    
    dst = pandas.DataFrame()
    dst["Q"] = dataset["Q"].copy()
    dst["I"] = np.full(dataset["Q"].shape, np.nan) # initialize as nan
    dst["IError"] = np.full(dataset["Q"].shape, np.nan)
    dst["Mask"] = np.zeros(dataset["Q"].shape, dtype = bool) # none masked

    if mirror:
        dataset = (dataset.iloc[::-1]).append(dataset, ignore_index = True)
        dataset.Q[0:dst.Q.size] *= -1
    if modType == "multiply":
        newQ = dataset["Q"] * QModifier
    elif modType == "shift":
        newQ = dataset["Q"] + QModifier
    elif modType == "slitSmear":
        newQ = np.sign(dataset.Q) * np.sqrt(dataset.Q**2 - QModifier**2)
    elif modType == "invSlitSmear":
        newQ = np.sign(dataset.Q) * np.sqrt(dataset.Q**2 + QModifier**2)
    else:
        print("modType not understood, should be 'multiply' or 'shift'.")

    
    # interpolator (linear) to equalize Q. 
    fI = interp1d(newQ, dataset["I"],
            kind = "linear", bounds_error = False)
    fE = interp1d(newQ, dataset["IError"],
            kind = "linear", bounds_error = False)

    # interpolate, rely on Mask to deliver final limits
    dst["I"] = fI(dst["Q"])
    dst["IError"] = fE(dst["Q"])

    # extra mask clip based on I or IError values:
    dst["Mask"] |= dst["I"].isnull()
    dst["Mask"] |= dst["IError"].isnull()
    return dst


