# -*- coding: utf-8 -*-

import vtk

def getSTLReader(filename): 
    # sets up a read function with the STL filename in place. 
    # Get an STL reader (https://pyscience.wordpress.com/2014/09/21/ray-casting-with-python-and-vtk-intersecting-linesrays-with-surface-meshes/):
    readerSTL = vtk.vtkSTLReader()
    readerSTL.SetFileName(filename)
    # 'update' the reader i.e. read the .stl file
    readerSTL.Update()
    return readerSTL

def STLToPolydata(filename): 
    # converts the STL data to polydata, which is what we need. 
    readerSTL = getSTLReader(filename)
    polydata = readerSTL.GetOutput()

    # If there are no points in 'vtkPolyData' something went wrong
    if polydata.GetNumberOfPoints() == 0:
        raise ValueError(
            "No point data could be loaded from '" + filename)
        return None
    # be nice:
    polydata.ComputeBounds()
    return polydata
