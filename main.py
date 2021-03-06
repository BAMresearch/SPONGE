# copied from the main.py file of the old McSAS:

import argparse
import logging
import multiprocessing
import sys, os
from sys import platform

# adapted from: https://stackoverflow.com/questions/8220108/how-do-i-check-the-operating-system-in-python
def isLinux():
    return platform == "linux" or platform == "linux2"
def isMac():
    return platform == "darwin"
def isWindows():
    return platform == "win32"

# from log import replaceStdOutErr
# from utils import isMac, isLinux

def getScriptPath():
    """Returns the full path to the current script file which calls this
    function."""
    thisFile = sys.executable # works for frozen app
    try: # __file__ not defined if frozen
        thisFile = os.path.join(os.getcwd(), __file__)
    except NameError: pass
    thisFile = os.path.abspath(thisFile)
    path, fn = os.path.split(thisFile)
    if os.path.isfile(path):
        path = os.path.dirname(path)
    return path, fn

# get script/executable location, add to module search path
SCRIPT_PATH, SCRIPT_FILENAME = getScriptPath()
if not hasattr(sys, "frozen") and SCRIPT_PATH not in sys.path:
    sys.path.append(SCRIPT_PATH)
    # FIXME? sphinx assumes the upper directory as root,
    # prefixes every module path with *mcsas*
    # -> do it the same way? Would simplify some things ...

if (isLinux() and hasattr(sys, "frozen")
    and "LD_LIBRARY_PATH" not in os.environ):
    os.environ["LD_LIBRARY_PATH"] = SCRIPT_PATH

import s


def makeAbsolutePath(relpath):
    return os.path.abspath(os.path.join(SCRIPT_PATH, relpath))

def main(argv = None):
    parser = argparse.ArgumentParser(description = """
            Simulates small-angle scattering patterns from 
            STL-descriptions of object surfaces. Can include
            polydispersity in size (uniformly scaling in all dimensions)
            Cobbled together by Brian R. Pauw.
            Released under a GPLv3+ license.
            """)
    # TODO: add info about output files to be created ...
    parser.add_argument("-f", "--efName", type = str, default = None,
            help = "Path to excel filename containing the sim settings")
    parser.add_argument("-g", "--group", type = str, default = None,
            help = "simulation group to work on")    
    parser.add_argument("-n", "--numProcesses", type = int, default = None,
            help = "Maximum number of parallel processes. Default is number of cores on system.")    
    parser.add_argument("-s", "--filename", type = str, default = None,
            help = "FOR DIRECT CALL: input STL filename")    
    parser.add_argument("-q", "--qmin", type = float, default = 0.01,
            help = "FOR DIRECT CALL: minimum q value")    
    parser.add_argument("-Q", "--qmax", type = float, default = 2,
            help = "FOR DIRECT CALL: maximum q value")    
    parser.add_argument("-N", "--nq", type = int, default = 100,
            help = "FOR DIRECT CALL: number of q values")    
    parser.add_argument("-R", "--nrep", type = int, default = 100,
            help = "FOR DIRECT CALL: number of repetitions")
    parser.add_argument("-P", "--npoints", type = int, default = 1000,
            help = "FOR DIRECT CALL: number of points in object")
    parser.add_argument("-M", "--mu", type = float, default = 1,
            help = "FOR DIRECT CALL: mean of Gaussian size distribution")
    parser.add_argument("-S", "--sigma", type = float, default = 0.01,
            help = "FOR DIRECT CALL: width (sigma) of Gaussian size distribution")
    parser.add_argument("-O", "--ofname", type = str, default = None,
            help = "FOR DIRECT CALL: output filename")
    parser.add_argument("-D", "--projectDirectory", type = str, default = '.',
            help = "FOR DIRECT CALL: project directory (starting point)")
 
    if isMac():
        # on OSX remove automatically provided PID,
        # otherwise argparse exits and the bundle start fails silently
        for i in range(len(sys.argv)):
            if sys.argv[i].startswith("-psn"): # PID provided by osx
                del sys.argv[i]
    try:
        args = parser.parse_args()
    except SystemExit as e:
        # useful debugging code, ensure destination is writable!
        # logfn = ("/tmp/{name}_unsupported_args.log"
        #   .format(name = SCRIPT_FILENAME))
        # with open(logfn, "w") as fd:
        #   fd.write("argv: " + str(sys.argv) + "\n")
        raise
    # initiate logging (to console stderr for now)
    # replaceStdOutErr() # replace all text output with our sinks
    adict = vars(args)
    s.s(**adict)

if __name__ == "__main__":
    multiprocessing.freeze_support()
    #manager=pyplot.get_current_fig_manager()
    #print manager
    #process input arguments
    main()
