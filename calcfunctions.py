# -*- coding: utf-8 -*-

# use in an instance with VTK!

import vtk 
import numpy as np
import scipy.spatial
import logging

def polydataToMass(polydata):
    # returns the object volume
    Mass = vtk.vtkMassProperties()
    Mass.SetInputData(polydata)
    Mass.Update() 
    return Mass.GetVolume()

def polydataToSurface(polydata):
    # returns the surface area
    Mass = vtk.vtkMassProperties()
    Mass.SetInputData(polydata)
    Mass.Update() 
    return Mass.GetSurfaceArea()

def pickPointsInMeshV2(mesh, nPoints = 1000): 
    # choose points within the boundaries. These points are then (double)checked whether 
    # they lie inside or outside the object. Points outside the object are discarded. 
    # this process is repeated until nPoints have been found inside. 
    
    # Find the limits of the mesh:
    mesh.ComputeBounds() # already done upon STL read
    (xMin, xMax, yMin, yMax, zMin, zMax) = mesh.GetBounds()
    # print("Limits: x: {}, {}, y: {}, {}, z: {}, {}".format(xMin, xMax, yMin, yMax, zMin, zMax))

    nFound = 0
    pts = []
    TPCoord = np.zeros([nPoints, 3]) # test block
    # inCoord = np.zeros([nPoints, 3]) # final point set
    while (nFound < nPoints):
        # generate a block of points to test
        TPCoord[:, 0] = np.random.uniform(low = xMin, high = xMax, size = nPoints)
        TPCoord[:, 1] = np.random.uniform(low = yMin, high = yMax, size = nPoints)
        TPCoord[:, 2] = np.random.uniform(low = zMin, high = zMax, size = nPoints)
        # add to vPts object:
        # version using vtk points object:

        TPts = vtk.vtkPoints()
        TPts.SetDataType(vtk.VTK_DOUBLE)
        dummy = [TPts.InsertNextPoint([TPCoord[j, 0], TPCoord[j, 1], TPCoord[j, 2]]) for j in range(nPoints)]
        chkPts = vtk.vtkPolyData()
        chkPts.SetPoints(TPts)

        # set up location checker, parts of this may be moved outside loop later:
        sel = vtk.vtkSelectEnclosedPoints()
        sel.SetInputData(chkPts)
        sel.SetSurfaceData(mesh)
        sel.CheckSurfaceOn()
        sel.Update()
        
        pointi = [] # new list
        j = 0
        while (nFound < nPoints) and (j < nPoints):
            if sel.IsInside(j):
                pointi.append(j)
                nFound += 1
            j+=1
            
        # add to final set:
        [pts.append(chkPts.GetPoint(j)) for j in pointi]
        # for j in pointi:
        #     inCoord[j,:] = chkPts.GetPoint(j)
        # print("{} points found of requested {}".format(nFound, nPoints))

    return pts

def logEdges(dist, qmin, qmax, nq):
    """
    Calculates the optimal histogramming bin edges for pointsToScatterD, based on input arguments:
    * dist *: the point-to-point distance list from scipy.spatial.distance.pdists
    * qmin *: the minimum requested Q value
    * qmax *: the maximum requested Q value
    * nq *: the number of requested q points 
    """
    logging.error('This "fast method"-functionality is depreciated, please use Debyer for this approach')

    # # core logEdges:    
    # logEdges = np.logspace(
    #     np.log10(2 * np.pi / qmax), 
    #     np.log10(2 * np.pi / qmin), 
    #     nq)

    # # if there are distances below:
    # if (dist.min() < (2 * np.pi / qmax)):
    #     logEdges = np.concatenate(
    #         [np.logspace(np.log10(dist.min()),np.log10(2 * np.pi / qmax),10, endpoint = False), 
    #         logEdges]
    #     )
    # # if there are distances above:
    # if (dist.max() > (2 * np.pi / qmin)):
    #     logEdges = np.concatenate(
    #         [logEdges, np.logspace(np.log10(2 * np.pi / qmin), np.log10(dist.max()),11, endpoint = True)[1:] 
    #         ]
    #     )
    # return logEdges

def pointsToScatterD(q, points, memSave = False):
    """ Calculate the scattering intensity from an array of points, by histogramming first. 
    This should be the quicker -- but potentially riskier -- method of calculating the 
    scattering intensity compared to the original pointsToScatter function. """

    logging.error('This "fast method"-functionality is depreciated, please use Debyer for this approach')
    # points = np.array(points)
    # dist = scipy.spatial.distance.pdist(points, metric = "euclidean")
    # lEdges = logEdges(dist, q.min(), q.max(), q.size)
    # dlog, elog = np.histogram(dist, bins = lEdges, density = False)
    # de = np.diff(elog) / 2 + elog[:-1] # middle distance of a given bin
    # # dle = dlog / np.diff(elog) # normalized fraction of contributions in a given bin
    # I2 = np.empty(q.shape)
    # I2.fill(np.nan) # initialize as nan
    # for qi, qval in enumerate(q):
    #     I2[qi] = 4 * np.pi * (dlog * np.sinc(de * qval / np.pi)).sum() / points.size**2
    # return I2

def pointsToScatter(q, points, memSave = False):
    # calculate the distance matrix between points, using a fast scipy function. 
    # This scipy function returns only unique distances, so only one distance 
    # value is returned for point1-point2 and point2-point1 combinations. It also
    # removes the zero distances between point1-point1. 
    # we then calculate the scattering using the Debye equation. 
    points = np.array(points)
    dist = scipy.spatial.distance.pdist(points, metric = "euclidean")
    if not memSave:
        inter = np.outer(np.abs(dist), q)
        # definition of np.sinc contains an additional factor pi, so we divide by pi. 
        # I = 2 * (np.sinc(inter / np.pi)).sum(axis=0) / points.size**2
        # prefactor should be 4 \pi.. perhaps. -> let's not for now. makes abs intensity scaling easier later
        I = (np.sinc(inter / np.pi)).sum(axis=0) / points.size**2
    else:
        I = np.empty(q.shape)
        I.fill(np.nan) # initialize as nan
        for qi, qval in enumerate(q):
            I[qi] = (np.sinc(dist * qval / np.pi)).sum() / points.size**2

    return I # , dist

