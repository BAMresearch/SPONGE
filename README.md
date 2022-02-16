# SPONGE
Simulates X-ray and Neutron scattering patterns from arbitrary shapes defined by STL files

# USAGE
## Command-line single simulation
The command-line is the recommended interface for using the SPONGE, as it will correctly exploit multiple threads. There are several command-line options, for example: 

'''
!python /Users/brian/Code/sponge/main.py \
-D /Users/brian/Code/sponge/ \
--filename "test/models/cubes/CubeL2.stl" \
-O "testDir/test.out" \
-n 4 -N 400 -q 0.005 -Q 5 -R 25 -P 2000 -M 2 -S 0.02
'''

This runs the main code, and starts from the working directory /Users/brian/Code/sponge/. 
The input STL file resides in the subdirectory (from the working directory) of tests/models/cubes/CubeL2.stl, a cube with a side of two nm. 
The result will be stored in the output file as three-column ascii, alongside a NeXus/HDF5 file with the simulation details. 
This will run on four cores, calculates 400 points in Q-space, ranging from 0.005 to 5 1/nm. 
It will calculate 25 independent repetitions, with 2000 points inside the object surface. 
The object will be scaled with scaling factors determined by a Gaussian distributino with mean mu=2 and sigma=0.02. 

## multiple simulations defined by an excel-sheet of parameters
[description will be added later]

## run from a Jupyter notebook
Note: when run from inside a Jupyter notebook, this will run on a single core only!
[description will be added later]
