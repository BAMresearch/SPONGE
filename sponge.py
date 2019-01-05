# -*- coding: utf-8 -*-

import pandas
import numpy as np

from .plotfunctions  import simPlot
from .calcfunctions  import pickPointsInMeshV2, pointsToScatter, pointsToScatterD, logEdges, polydataToMass
from .stlfunctions   import getSTLReader, STLToPolydata
from .distfunctions  import interpolate, distSpreadGaussian
from .smearfunctions import halfTrapzPDF, slitSmearTrapz

from multiprocessing import Pool

class sponge(object):
    # generates an instance of Sponge
    # initialise
    # q = None            # q vector
    # mesh = None         # the STL mesh
    # volume = None       # object volume, used for final scaling
    rawData = list()    # list of data results of individual reps
    data = dict()       # in this sequence: averaged over runs
    distData = dict()   # smeared for dilation, size distribution
    smearData = dict()  # smeared for beam shape (slit only)
    parameters = {      # some calculation parameters
            "filename"     : "SphereR1F200.stl",  # input shape, STL
            "qMin"         : 0.1,                 # min. q to calc
            "qMax"         : 20,                  # max. q to calc
            "nRep"         : 10,                  # independent repetitions
            "nq"           : 100,                 # number of q points
            "nPoints"      : 1000,                # points in object
            "smear"        : False,               # slit smearing?
            "dist"         : False,               # size distribution?
            "sigma"        : 0.2,                 # size dist. width (Gaussian)
            "ndiv"         : 20,                  # number of integration divs
            "halfUmbra"    : 1.5,                 # slit smear width
            "halfPenUmbra" : 2.5,                 # slit smear width
            "parallel"     : False,               # try running this in parallel
            "threads"      : 3,                   # with this many threads
            "fastMethod"   : True                 # "fastMethod" uses pointsToScatterD instead of pointsToScatter
            } 

    def __init__(self, params = None):
        # reset everything for a new instance
        self.rawData = list()    # list of data results of individual reps
        self.data = dict()       # in this sequence: averaged over runs
        self.distData = dict()   # smeared for dilation, size distribution
        self.smearData = dict()  # smeared for beam shape (slit only)
        if params is not None:
            self.parameters.update(params)
        # initialize Q:
        self.q = np.logspace(
                np.log10(self.parameters["qMin"]), 
                np.log10(self.parameters["qMax"]), 
                self.parameters["nq"])

    # replaced in notebook..
    
    def singleRun(self, position = None):
        """ Starts a single calculation based on the provided parameters """
        mesh = STLToPolydata(self.parameters["filename"])
        pts = pickPointsInMeshV2(mesh, self.parameters["nPoints"])
        if not self.parameters["fastMethod"]:
            I = pointsToScatter(self.q, pts)
        else:
            I = pointsToScatterD(self.q, pts)
        self.rawData.append(I)
        return 
        
    def multiRun(self):
        # repeats the intensity calculation a number of times to get a good average intensity, 
        # and get uncertainty estimates to boot. 
        self.rawData = list()
        
        if self.parameters["parallel"]: # doesn't work yet...
            with Pool(self.parameters["threads"]) as p:
                result = p.map(self.singleRun, [range(self.parameters["nRep"])])
                p.close()
                p.join()
        else:
            for rep in range(self.parameters["nRep"]):
                self.singleRun()

        self.data = pandas.DataFrame({"Q": self.q})
        self.data["I"] = np.array(self.rawData).mean(axis = 0)
        self.data["IError"] = np.array(self.rawData).std(axis = 0, ddof = 1) / np.sqrt(np.array(self.rawData).shape[0])
        return

    def applySizeDistribution(self):
        x, g, addDat = distSpreadGaussian(
                self.data, 
                sigma = self.parameters["sigma"], 
                ndiv = self.parameters["ndiv"])
        addDat.update({
                "distX" : x,
                "distPX" : g})
        self.distData = pandas.DataFrame(addDat)

        return 

    def applySmearing(self):
        x, g, addDat = slitSmearTrapz(
                self.distData, 
                halfUmbra = self.parameters["halfUmbra"], 
                halfPenUmbra = self.parameters["halfPenUmbra"],
                ndiv = self.parameters["ndiv"])
        addDat.update({
            "smearX" : x,
            "smearPX" : g})
        self.smearData = pandas.DataFrame(addDat)

        return 


