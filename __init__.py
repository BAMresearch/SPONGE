# -*- coding: utf-8 -*-
# __init__.py

"""
General imports of sponge
"""

__author__ = "Brian R. Pauw"
__contact__ = "brian@stack.nl"
__license__ = "GPLv3+"
__copyright__ = "Bundesanstalt für Materialforschung und -Prüfung"
__date__ = "2019-06-21"
__status__ = "beta"
version = "17"

__all__ = []

from .sponge import sponge
from .plotfunctions import simPlot
from .calcfunctions import pickPointsInMeshV2, pointsToScatter, pointsToScatterD, logEdges, polydataToMass
from .stlfunctions  import getSTLReader, STLToPolydata
from .distfunctions import interpolate, distSpreadGaussian
from .smearfunctions import halfTrapzPDF, slitSmearTrapz

def argparser():
    parser = argparse.ArgumentParser(description = """
            Simulates small-angle scattering patterns from 
            STL-descriptions of object surfaces. Can include
            polydispersity in size (uniformly scaling in all dimensions)
            Cobbled together by Brian R. Pauw.
            Released under a GPLv3+ license.
            """)
    parser.add_argument("-f", "--efName", type = str, default = None,
            help = "Path to excel filename containing the sim settings")
    parser.add_argument("-g", "--group", type = str, default = None,
            help = "simulation group to work on")    
    return parser.parse_args()

if __name__ == "__main__":
    #manager=pyplot.get_current_fig_manager()
    #print manager
    #process input arguments
    adict = argparser()
    #run the program, scotty! I want a kwargs object, so convert args:
    adict = vars(adict)
    sponge.sponge(**adict) #and expand to kwargs

def "__main__":
    print("boo")


# vim: set ts=4 sts=4 sw=4 tw=0:
