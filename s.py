# -*- coding: utf-8 -*-

import pandas
import numpy as np
import argparse
import datetime
import os

from .nexuswriter import NeXusWriter

from pathlib import Path

# from plotfunctions  import simPlot
from . import calcfunctions
# from calcfunctions  import pickPointsInMeshV2, pointsToScatter, pointsToScatterD, logEdges, polydataToMass
from . import stlfunctions
# from stlfunctions   import getSTLReader, STLToPolydata
# from .distfunctions  import interpolate, distSpreadGaussian
# from .smearfunctions import halfTrapzPDF, slitSmearTrapz
import multiprocessing 
# from multiprocessing import Pool

class s(object):
    # generates an instance of Sponge
    resultDict = None   # results in dictionary form
    volumes = None      # volumes of the simulated objects
    surfaceAreas = None # surface areas of the simulated objects
    simdata = None      # Q, I, IError
    IAtQMin = None      # I at smallest q, not volume-square compensated, should be close to 1 or you didn't simulate to low enough q
    _SLDs = None
    _volumefraction = None

    def __init__(self, numProcesses = 1, run=True, **kwargs):
        """ 
        can be run either using an excel file as input or for direct simulation of a single file.
        
        
        """
        # reset everything for a new instance
        self.resultDict = []   # results in dictionary form
        self.volumes = []      # volumes of the simulated objects
        self.surfaceAreas = [] # surface areas of the simulated objects
        self.simdata = []      # Q, I, IError
        self.IAtQMin = []      # I at smallest q, not volume-square compensated, should be close to 1 or you didn't simulate to low enough q
        # array of SLDs for each of the phases, also for multiphase modelling: 
        # in the order of [SLD(phase1), SLD(phase2), ..., SLD(dispersant)]). 
        # For single-phase, just have [SLD(particle), SLD(dispersant)]
        self._SLDs = [1., 0.]
        self._volumefraction = 1. # volume fraction of scatterers        
        
        assert (kwargs['efName'] is not None) or (kwargs['filename'] is not None), 'either an excel file with a list of simulations or a single STL file should be provided'
        if kwargs['efName'] is not None: 
            # run in the traditional way:
            if run:
                assert kwargs['efName'] is not None, 'Path to excel filename with the simulation settings must be provided'
                assert kwargs['group'] is not None, 'simulation group must be specified'
                Tests = self.loadTests(kwargs['efName'])
                self.resultDict = self.runTests(Tests = Tests, group = kwargs['group'], numProcesses = numProcesses)
        else: # we run this in the single-file way:
            if run:
                Tests = self.genTest(kwargs)
                resultDict = self.runTests(Tests=Tests, group = 'default',numProcesses = numProcesses)
            
    def loadTests(self, efName):
        efName = Path(efName) # it's ok if this is done multiple times...
        Tests = pandas.DataFrame() # CLEARS TESTS!

        # efName = Path('Martin/Martin.xlsx')
        projectdirectory = efName.parent
        excelFilename = efName.name

        Tests = pandas.read_excel(efName, skiprows = 1)
        Tests = Tests.dropna(axis = 0, how = "all") #remove empty rows for cleaning up. 
        Tests.columns = Tests.columns.str.lower() # lower case column names only
        # os.chdir(kansas)
        # cast to the right datatypes:
        Tests = Tests.astype({
            "npoints":"int", 
            "nq":"int", 
            "nrep":"int", 
            "memsave":"bool", 
            # "fastmethod":"bool",
            "qmin": "float",
            "qmax": "float",
            "mu": "float",
            "sigma": "float",
            })
        # Tests['fastmethod'] = False # no longer implemented...  
        Tests["projectdirectory"] = projectdirectory # add a project directory to all the entries  
        return Tests

    def genTest(self, pars):
        """Similar to loadTests, but this is generating test parameters for a single simulation based on input parameters"""
        # fill defaults:
        Tests = pandas.DataFrame(
            data={
                "testgroup": "default",
                "filename": Path('.'),
                "projectdirectory": Path('.'),
                "ofname": Path('default.out'),
                "npoints":1000, 
                "nq":200, 
                "nrep":100, 
                "memsave":True, 
                # "fastmethod":False, # no longer implemented
                "qmin": 0.01,
                "qmax": 2,
                "mu": 1,
                "sigma": 0.01,
            }, 
            index = [0],
        )
        for key, val in pars.items():
            if key.lower() in ['filename', 'ofname', 'projectdirectory']:
                val = Path(val)
            print(f'Changing parameter: {key.lower()} to value: {val}')
            Tests[key.lower()] = val
        # cast to the right datatypes:
        Tests = Tests.astype({
            # "filename": Path,
            # "ofname": Path,
            # "projectdirectory": Path,
            "npoints":"int", 
            "nq":"int", 
            "nrep":"int", 
            "memsave":"bool", 
            "qmin": "float",
            "qmax": "float",
            "mu": "float",
            "sigma": "float",
            })
        return Tests

    def singleRun(self, parameters):
        """ Starts a single calculation based on the provided parameters """
        q = np.logspace(
            np.log10(parameters["qmin"]),
            np.log10(parameters["qmax"]),
            parameters["nq"])

        mesh = stlfunctions.STLToPolydata(Path(parameters["projectdirectory"], parameters["filename"]).as_posix())
            
        pts = calcfunctions.pickPointsInMeshV2(mesh, parameters["npoints"])
        # make sure we don't get too negative. 
        scaler = -1.0
        while (scaler < 0):
            scaler = np.random.normal(loc = parameters["mu"], scale = parameters["sigma"]) # scaling factor to apply to the shape
        
        # if not parameters["fastmethod"]: # fastMethod depreciated
        # this is the slowest I think so far. 6s with memmsave, 25s without (memory IO problem)
        I = calcfunctions.pointsToScatter(q * scaler, pts, parameters["memsave"])
        # else:
        #     I = calcfunctions.pointsToScatterD(q * scaler, pts, parameters["memsave"])

        vol = calcfunctions.polydataToMass(mesh) * scaler**3 # I think this is correct with the scaler
        # bonus surface area (for surface-to-volume ratios):
        surf = calcfunctions.polydataToSurface(mesh) * scaler **2

        return I, vol, surf

    def multiRun(self, parameters, numProcesses = None):
        if "sld" in parameters:
            self._SLD = list(parameters["sld"])
        if "volumefraction" in parameters:
            self._volumefraction = float(parameters["volumefraction"])

        q = np.logspace(
            np.log10(parameters["qmin"]),
            np.log10(parameters["qmax"]),
            parameters["nq"])

        if numProcesses is None:
            nump = multiprocessing.cpu_count()
        else:
            assert isinstance(numProcesses, int)
            nump = np.minimum(multiprocessing.cpu_count(), numProcesses)
        if nump > 1:
            # multiprocessing
            Pool = multiprocessing.Pool(processes = nump)
            mapParam = [parameters for i in range(int(parameters["nrep"]))]
            rawData = Pool.map(self.singleRun, mapParam)    
            Pool.close()
            Pool.join()
        else:
            # single threaded:
            rawData = []
            for _ in range(int(parameters["nrep"])):
                rawData.append(self.singleRun(parameters))
        
        # pick apart intensities and volume outputs:
        rawDataI = []
        rawIAtQMin = [] # this should be close to 1, otherwise the sim wasn't done to low enough q to ge ta good Guinier
        rawDataV = []
        rawDataS = [] # surface area
        for item in rawData:
            # note: we are volume-weighting the intensity here!
            # if you want "standard" number-weighted intensity, then multiply not by item[1] but item[1]**2
            rawDataI.append(item[0] * item[1] * self._volumefraction * (self._SLDs[-1] - self._SLDs[0])**2) # =Fsq/V of SasView-style (intensity multiplied by volume-square, then divided by volume once). Still should be multiplied by deltaSLD**2 to go to abs units, scale will be volume fraction if obtained intensity from the calculation converges to 1 at q=0
            rawIAtQMin.append(item[0][0])
            # we have surfaces available at item[2] if we want. 
            rawDataV.append(item[1]) # volumes that come out are actual volumes (e.g. for sphere = 4/3 pi r^3)
            rawDataS.append(item[2]) # volumes that come out are actual volumes (e.g. for sphere = 4/3 pi r^3)
        rawDataI = np.array(rawDataI)
        self.IAtQMin = np.array(rawIAtQMin)
        self.volumes = np.array(rawDataV)
        self.surfaceAreas = np.array(rawDataS)

        data = pandas.DataFrame({"Q": q}) 
        data["I"] = rawDataI.mean(axis = 0) # * rawDataV.mean() # second half of scaling, first half is done in I.
        data["IError"] = rawDataI.std(axis = 0, ddof = 1) / np.sqrt(rawDataI.shape[0]) # * rawDataV.mean()
        # also need to save volumes and surfaces somewhere. 
        self.simdata = data

        if parameters["ofname"] is not None:
            Path(parameters["projectdirectory"],parameters["ofname"]).parent.mkdir(parents=True, exist_ok=True)
            data.to_csv(Path(parameters["projectdirectory"],parameters["ofname"]).as_posix(), header = False, sep = ';', index = False)

        return {"data"      : data,
               "parameters" : parameters}

    def runTests(self, Tests = None, start = 0, stop = None, group = None, numProcesses = None):
        """Runs a series of simulations as defined in the Tests loaded from the excel files"""
        resultDict = {}
        print(Tests)
        if stop is None:
            testindices = Tests.index.values
        if group is not None:
            testindices = Tests[Tests.testgroup == group].index.tolist() # old: .values
                    
        for tn, testindex in enumerate(testindices):
            print("Test: {} of {}".format(tn + 1, len(testindices)))
            param = Tests.loc[testindex]
            try:
                del res
            except NameError:
                pass
            except:
                raise

            res = self.multiRun(param, numProcesses = numProcesses)
            self.storeResult(res["parameters"], res["data"])
            resultDict.update({testindex: res})
        return resultDict

    def storeResult(self, parameters, data):
        """Stores the result in a NeXus structure... """
        directDict = {}
        for rdKey, rdValue in parameters.items():
            if rdKey == 'data':
                continue # skip this, not part of the metadata
            directDict.update(
                {'/sasentry1/simulationParameters/{}'.format(rdKey) : '{}'.format(rdValue)}        
            )
        for rdKey in [
            'volumes',
            'surfaceAreas',
            'simdata',
            'IAtQMin',
            ]:
            directDict.update(
                {f'/sasentry1/simulationMetaValues/{rdKey}' : f'{getattr(self, rdKey)}'}        
            )
        
        tp = Path(parameters['projectdirectory'], parameters['ofname'])
        NeXusWriter(
            filename = tp.with_suffix('.nxs'), 
            Q = data['Q'], 
            I = data['I'],
            IError = data['IError'], 
            wavelength = 0.1542,
            wavelengthUnits = 'nm',
            title = 'Simulated data for testgroup: {}, file: {}'.format(
                parameters['testgroup'], parameters['filename']            
            ),
            timestamp = datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc).replace(microsecond=0).isoformat(),
            overwrite = True,
            directDict = directDict, 
        )        

# To move the resulting intensity to absolute units, do the following:
# Multiply the intensity by the square of the SLD contrast between object and matrix
# multiply this with the volume fraction of the object in solution