import sys

sys.path.append("./pcrasterModules/") # from PCRasterModules - still needed for other scripts to work?

from pcrasterModules import generalfunctions_test01
from pcraster.framework import *

aVariable = generalfunctions_test01.openSamplesAndTimestepsAsNumpyArraysAsNumpyArray(
    'biTS', 30, range(100, 5200, 100))
numpy.save('biTS', aVariable)
