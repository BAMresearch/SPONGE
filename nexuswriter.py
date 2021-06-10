# -*- coding: utf-8 -*-
# nexuswriter.py
# Author: Brian R. Pauw
# License: GPLv3

# helper function for writing datasets:
import os
import h5py
import numpy as np

class NeXusWriter(object):
    """ 
    Writes the provided data and keys into a NeXus-/NXcanSAS-conform structure in an HDF5 file
    """
    
    _attrDict = {}
    _filename = None
    _overwrite = None
    _inputMapping = {}
    _requiredInputKeys = ["Q", "I", "IError"]

    ######## INPUT KEYWORD-VALUE MAPPING IS DEFINED HERE #########
    """ 
    The following dictionary identifies possible input keywords to NeXusWriter, 
    and their respective mapping paths to datasets. 
    Any mapping paths containing an '@'-symbol will link to the attribute table. In 
    that case, the element of the mapping path after the '@'-symbol is assumed to be the 
    name of the attribute
    """
    _inputMapping = {
        "Q"          : "/sasentry1/sasdata1/Q", # dataset
        "I"          : "/sasentry1/sasdata1/I", # dataset
        "IError"     : "/sasentry1/sasdata1/Idev", # dataset
        "Qdev"       : "/sasentry1/sasdata1/Qdev", # optional dataset
        "wavelength" : "/sasentry1/instrument/incident_wavelength",
        "wavelengthUnits" : "/sasentry1/instrument/incident_wavelength@units",
        "sampleName" : "/sasentry1/sample/name", # dataset
        "title"      : "/sasentry1/title", # dataset
        "instrumentTitle" : "/sasentry1/instrument/title", # dataset
        "timestamp"  : "/sasentry1/sasdata1@timestamp", # attribute
        "QUnits"     : "/sasentry1/sasdata1/Q@units", # attribute
        "IUnits"     : "/sasentry1/sasdata1/I@units",
        "IdevUnits"  : "/sasentry1/sasdata1/Idev@units",
        "QdevUnits"  : "/sasentry1/sasdata1/Qdev@units",
    }

    ######## ATTRIBUTES ARE DEFINED HERE #########

    """Fills the attribute dictionary with fields we need"""
    _defaults = {
        # /
        "/@canSAS_class": "SASroot",
        "/@default": "sasentry1", 
        # /sasentry1
        "/sasentry1@NX_class": "NXentry", 
        "/sasentry1@canSAS_class": "SASentry", 
        "/sasentry1@version": "1.0", 
        # /sasentry/sasdata1
        "/sasentry1/sasdata1@NX_class": "NXdata", 
        "/sasentry1/sasdata1@canSAS_class": "SASdata", 
        "/sasentry1/sasdata1@I_axes": "Q", 
        "/sasentry1/sasdata1@Q_indices": "0", 
        "/sasentry1/sasdata1@timestamp": "2019-01-30T16:19:54+00:00", 
        "/sasentry1/sasdata1@signal": "I", 
        # /sasentry1/sasdata1/I
        "/sasentry1/sasdata1/I@units": "1/m/sr", 
        "/sasentry1/sasdata1/I@uncertainties": "Idev", 
        # /sasentry1/sasdata1/Idev
        "/sasentry1/sasdata1/Idev@units": "1/m/sr", 
        # /sasentry1/sasdata1/Q
        "/sasentry1/sasdata1/Q@units": "1/nm", 
        # /sasentry1/instrument
        "/sasentry1/instrument@canSAS_class": "SASinstrument",
        "/sasentry1/instrument@NX_class": "NXdata", 
        # /sasentry1/instrument/incident_wavelength
        "/sasentry1/instrument/incident_wavelength@units": "nm", 
        # /sasentry1/sample
        "/sasentry1/sample@canSAS_class": "SASsample", 
        "/sasentry1/sample@NX_class": "NXsample", 
    }


    def __init__(self, filename = None, overwrite = False, directDict = {}, **kwargs):
        """
        NeXusWriter is a flexible class for constructing, or adding to (HDF5-based) 
        NeXus files. It can be used using (a combination of) two mechanisms:
          1. A mapping can be defined in the NeXusWriter.inputMapping dictionary, which
             defines what additional input arguments should end up where in the NeXus file. 
             This can be convenient when using the NeXusWriter in other code which always
             writes the same standard output or should be fed with user-intelligible 
             parameters, 
             
             and/or
             
          2. A dictionary ("directDict") can be provided to the arguments of NeXusWriter,
             with a pathKey (e.g. "/entry1/instrument/detector/data" for a dataset, or 
             "entry1/sample/transmission@units" for an attribute), and a value to associate
             with it. 
        After the input arguments are handled, the defaults in the NeXusWriter._defaults 
        dictionary are filled in for all still missing entries. 
        
        NeXusWriter can be called with the following arguments:
          * filename *: (required) a valid filename or Path object
          * overwrite *: will delete and recreate if the file exists
          * directDict *: A dictionary containing 'pathKey: value'-combinations, that adhere to the specifications in 
        the addAttribute and addDataset methods.
          further keyword-value pairs can be defined later on in this function for default mappings.
          See NeXusWriter._inputMapping and NeXusWriter._defaults dictionaries for the current defaults


        """
        self._filename = filename
        self._overwrite = overwrite
        
        for key in self._requiredInputKeys:
            assert key in kwargs.keys(), "{} must be provided".format(key)
        
        # ensure we're good to go
        self.validate()
        
        self._dealWithInput(kwargs)
        self._dealWithInput(directDict)
        
        # find the attributes which haven't been defined yet, and fill them with the defaults
        for key, value in self._defaults.items():    
            if (key in kwargs.keys()) or (key in directDict.keys()):
                continue
            else:
                # place the remaining default attributes in the file:
                self.addAttribute(pathKey = key, value = value)
           
    def _dealWithInput(self, inputDict):
        """
        deals with the input arguments to the class that need to be written to the NeXus structure. 
        inputDict must contain pathKey: value combinations, that adhere to the specifications in 
        the addAttribute and addDataset methods.
        """
        # datasets must be set first, then attributes can be attached:
        for datasetsOrAttributes in ['datasets', 'attributes']:
            for key, value in inputDict.items():
                if not key in self._inputMapping:
                    # print("input key {} not in input mapping, interpreting directly...")
                    kv = [key, value]
                else:
                    kv = [self._inputMapping[key], value]
                if ('@' in kv[0]) and (datasetsOrAttributes.lower() == 'attributes'): # we are dealing with an attribute here
                    # print('adding attribute {}'.format(kv[0]))
                    self.addAttribute(pathKey = kv[0], value = kv[1])
                if not('@' in kv[0]) and (datasetsOrAttributes.lower() == 'datasets'): # we are dealing with an attribute here
                    # print('adding dataset {}'.format(kv[0]))
                    self.addDataset(pathKey = kv[0], value = kv[1])
                
    def validate(self):
        """Generic validation method"""
        assert self._filename is not None, "Output filename (filename) must be provided!"
        if not self._overwrite:
            assert not os.path.isfile(self._filename), "Filename cannot exist already"
        else:
            if os.path.isfile(self._filename):
                os.remove(self._filename)
        # additional checks can be put here too...
        
    def addAttribute(self, pathKey = None, value = None):
        """
        Adds attributes to a dataset or group in NeXus. 
        pathKey must be corresponding to:
        /path/to/datasetOrGroup@attribute
        value can be any valid dtype.         
        """
        # adapted from Structurize:
        assert pathKey is not None, "HDF5 path/key combination cannot be empty"
        print("adding attribute at pathKey {}".format(pathKey))
        path, key = pathKey.rsplit("@", 1)
        
        with h5py.File(self._filename, 'a') as h5f:
            if not path in h5f:
                print("    Location {} does not exist in output file".format(path))
                return

            # write attribute
            h5f[path].attrs[key] = value
        
    def addDataset(self, pathKey = None, value = None):
        """
        Adds (or replaces) a dataset to a NeXus structure. 
        pathKey must be corresponding to:
        /path/to/dataset
        value can be any valid dtype. 
          - Lists are converted to np.arrays.
          - Unicode string arrays are converted to h5py special_dtype '<U6'
        """
        # adapted from McSAS3:
        assert pathKey is not None, "HDF5 path/key combination cannot be empty"
        print("adding dataset at pathKey {}".format(pathKey))
    
        path, key = pathKey.rsplit('/', 1)
        
        """stores the settings in an output file (HDF5)"""
        with h5py.File(self._filename, 'a') as h5f:
            h5g = h5f.require_group(path)

            # store arrays:
            # convert all compatible data types to arrays:
            if type(value) is tuple or type(value) is list:
                value = np.array(value)
            if value is not None and type(value) is np.ndarray:
                # HDF cannot store unicode string arrays, these need to be stored as a special type:
                if value.dtype == '<U6':
                    value = value.astype(h5py.special_dtype(vlen=str))
                # store the data in the prefiously defined group:
                h5g.require_dataset(key, data = value, shape = value.shape, dtype = value.dtype, compression = "gzip")

            # non-array values are stored here:
            elif value is not None:
                # try and see if the destination already exists.. This can be done by require_dataset, but that requires shape and dtype to be specified. This method doesn't:
                dset = h5g.get(key, None)
                if dset is None:
                    h5g.create_dataset(key, data = value)
                else:
                    dset[()] = value

