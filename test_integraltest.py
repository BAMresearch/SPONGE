import unittest
import pandas, os
from s import s 
# some tests to ensure the functionality persists. This bypasses main.py, and directly addresses s
import warnings
import numpy as np

class testSponge(unittest.TestCase):
    def testSponge(self):
        # new instance:
        spongeInstance = s(run=False)
        # let's set up a set of testing parameters:
        tests = pandas.DataFrame(data= {
            'testgroup':'unittest',
            'filename': 'test/models/spheres/SphereR1F200.stl',
            'qmin': 0.1,
            'qmax': 160,
            'nq': 400, 
            'nrep': 10, # normally 100
            'npoints': 1000, # for testing only, you can go up to 32k, memory permitting
            'mu': 1., # mean scaling factor
            'sigma': 0., # (gaussian distribution width), now vol-weighted!
            'memsave': True, # speeds things up a little
            'ofname': 'test/simdata/sphere_unittest.dat',
            'projectdirectory': os.getcwd(),
        }, index = [0])
        res = spongeInstance.multiRun(tests.loc[0], numProcesses = 1) # numProcesses=1 runs without multiprocessing
        spongeInstance.storeResult(res["parameters"], res["data"])
        np.testing.assert_allclose(spongeInstance.IAtQMin, 0.998, atol = 0.001)
        np.testing.assert_allclose(spongeInstance.volumes, 4.19, atol = 0.01)

        spongeInstance2 = s(run=False)
        # let's set up a set of testing parameters:
        tests = pandas.DataFrame(data= {
            'testgroup':'unittest',
            'filename': 'test/models/spheres/SphereR3F200.stl',
            'qmin': 0.1,
            'qmax': 160,
            'nq': 400, 
            'nrep': 10, # normally 100
            'npoints': 1000, # for testing only, you can go up to 32k, memory permitting
            'mu': 1., # mean scaling factor
            'sigma': 0., # (gaussian distribution width), now vol-weighted!
            'memsave': True, # speeds things up a little
            'ofname': 'test/simdata/sphereR3_unittest.dat',
            'projectdirectory': os.getcwd(),
        }, index = [0])
        res = spongeInstance2.multiRun(tests.loc[0], numProcesses = 2) # numProcesses=1 runs without multiprocessing
        spongeInstance2.storeResult(res["parameters"], res["data"])
        np.testing.assert_allclose(spongeInstance2.IAtQMin, 0.982, atol = 0.002)
        np.testing.assert_allclose(spongeInstance2.volumes, 113.1, atol = 0.1)

        spongeInstance3 = s(run=False)
        # let's set up a set of testing parameters:
        tests = pandas.DataFrame(data= {
            'testgroup':'unittest',
            'filename': 'test/models/spheres/SphereR0p5F200.stl',
            'qmin': 0.1,
            'qmax': 160,
            'nq': 400, 
            'nrep': 10, # normally 100
            'npoints': 1000, # for testing only, you can go up to 32k, memory permitting
            'mu': 1., # mean scaling factor
            'sigma': 0., # (gaussian distribution width), now vol-weighted!
            'memsave': True, # speeds things up a little
            'ofname': 'test/simdata/sphereR0p5_unittest.dat',
            'projectdirectory': os.getcwd(),
        }, index = [0])
        res = spongeInstance3.multiRun(tests.loc[0], numProcesses = 2) # numProcesses=1 runs without multiprocessing
        spongeInstance3.storeResult(res["parameters"], res["data"])
        np.testing.assert_allclose(spongeInstance3.IAtQMin, 0.999, atol = 0.001)
        np.testing.assert_allclose(spongeInstance3.volumes, 0.524, atol = 0.001)


if __name__ == "__main__":
    unittest.main()