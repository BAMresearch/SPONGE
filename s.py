# -*- coding: utf-8 -*-

import pandas
import numpy as np
import argparse
import datetime

from nexuswriter import NeXusWriter

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
        Tests = self.loadTests(efName)
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
        Tests["projectDirectory"] = projectDirectory # add a project directory to all the entries  
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
        
        tp = Path(parameters['projectDirectory'], parameters['ofname'])
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
