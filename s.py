# -*- coding: utf-8 -*-

import pandas
import numpy as np
import argparse

from pathlib import Path

# from plotfunctions  import simPlot
import calcfunctions
# from calcfunctions  import pickPointsInMeshV2, pointsToScatter, pointsToScatterD, logEdges, polydataToMass
import stlfunctions
# from stlfunctions   import getSTLReader, STLToPolydata
# from .distfunctions  import interpolate, distSpreadGaussian
# from .smearfunctions import halfTrapzPDF, slitSmearTrapz
import multiprocessing 
# from multiprocessing import Pool

class s(object):
    # generates an instance of Sponge

    def __init__(self, efName = None, group = None):
        # reset everything for a new instance
        if efName is None:
            error('Path to excel filename with the simulation settings must be provided')
        if group is None:
            error('simulation group must be specified')
        Tests = loadTests(efName)
        Tests["projectDirectory"] = projectDirectory # add a project directory to all the entries
        resultDict = self.runTests(Tests = Tests, group = group)
            
    def loadTests(self, efName):
        efName = Path(efName) # it's ok if this is done multiple times...
        Tests = pandas.DataFrame() # CLEARS TESTS!

        # efName = Path('Martin/Martin.xlsx')
        projectDirectory = efName.parent
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
            "fastmethod":"bool",
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

        mesh = stlfunctions.STLToPolydata(Path(parameters["projectDirectory"], parameters["filename"]).as_posix())
            
        pts = calcfunctions.pickPointsInMeshV2(mesh, parameters["npoints"])
        # make sure we don't get too negative. 
        scaler = -1.0
        while (scaler < 0):
            scaler = np.random.normal(loc = parameters["mu"], scale = parameters["sigma"]) # scaling factor to apply to the shape
        
        if not parameters["fastmethod"]:
            I = calcfunctions.pointsToScatter(q * scaler, pts, parameters["memsave"])
        else:
            I = calcfunctions.pointsToScatterD(q * scaler, pts, parameters["memsave"])

        vol = calcfunctions.polydataToMass(mesh) * scaler**3 # I think this is correct with the scaler
        # print("Correcting for volume: {} nm^3".format(vol))
        # return I * vol **2
        return I, vol

    def multiRun(self, parameters, progressBar = True):
        q = np.logspace(
            np.log10(parameters["qmin"]),
            np.log10(parameters["qmax"]),
            parameters["nq"])

        Pool = multiprocessing.Pool(processes = multiprocessing.cpu_count())
        mapParam = [parameters for i in range(int(parameters["nrep"]))]
        rawData = Pool.map(self.singleRun, mapParam)    
        Pool.close()
        Pool.join()
        
        # pick apart intensities and volume outputs:
        rawDataI = []
        rawDataV = []
        for item in rawData:
            rawDataI.append(item[0] * item[1]) # volume-weighting
            rawDataV.append(item[1])
        rawDataI = np.array(rawDataI)
        rawDataV = np.array(rawDataV)
            
        data = pandas.DataFrame({"Q": q}) 
        data["I"] = rawDataI.mean(axis = 0) * rawDataV.mean() # second half of scaling, first half is done in I.
        data["IError"] = rawDataI.std(axis = 0, ddof = 1) / np.sqrt(rawDataI.shape[0]) * rawDataV.mean()
        
        if parameters["ofname"] is not None:
            data.to_csv(Path(parameters["projectDirectory"],parameters["ofname"]).as_posix(), header = False, sep = ';', index = False)

        return {"data"      : data,
               "parameters" : parameters}

    def runTests(self, Tests = None, start = 0, stop = None, group = None):
        resultDict = {}
        print(Tests)
        if stop is None:
            testindices = Tests.index.values
        if group is not None:
            testindices = Tests[Tests.testgroup == group].index.tolist() # old: .values
                    
        for testindex in testindices:
            print("Testindex: {} of {}".format(testindex, len(testindices)))
            param = Tests.loc[testindex]
            try:
                del res
            except NameError:
                pass
            except:
                raise

            res = self.multiRun(param)
            resultDict.update({testindex: res})
        return resultDict

        
        
    # def multiRun(self):
    #     # repeats the intensity calculation a number of times to get a good average intensity, 
    #     # and get uncertainty estimates to boot. 
    #     self.rawData = list()
        
    #     if self.parameters["parallel"]: # doesn't work yet...
    #         with Pool(self.parameters["threads"]) as p:
    #             result = p.map(self.singleRun, [range(self.parameters["nRep"])])
    #             p.close()
    #             p.join()
    #     else:
    #         for rep in range(self.parameters["nRep"]):
    #             self.singleRun()

    #     self.data = pandas.DataFrame({"Q": self.q})
    #     self.data["I"] = np.array(self.rawData).mean(axis = 0)
    #     self.data["IError"] = np.array(self.rawData).std(axis = 0, ddof = 1) / np.sqrt(np.array(self.rawData).shape[0])
    #     return

    # def applySizeDistribution(self):
    #     x, g, addDat = distSpreadGaussian(
    #             self.data, 
    #             sigma = self.parameters["sigma"], 
    #             ndiv = self.parameters["ndiv"])
    #     addDat.update({
    #             "distX" : x,
    #             "distPX" : g})
    #     self.distData = pandas.DataFrame(addDat)

    #     return 

    # def applySmearing(self):
    #     x, g, addDat = slitSmearTrapz(
    #             self.distData, 
    #             halfUmbra = self.parameters["halfUmbra"], 
    #             halfPenUmbra = self.parameters["halfPenUmbra"],
    #             ndiv = self.parameters["ndiv"])
    #     addDat.update({
    #         "smearX" : x,
    #         "smearPX" : g})
    #     self.smearData = pandas.DataFrame(addDat)

    #     return 
