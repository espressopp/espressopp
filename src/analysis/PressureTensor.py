"""
************************************
**PressureTensor** - Analysis
************************************

This class computes the pressure tensor of the system.
It can be used as standalone class in python as well as
in combination with the integrator extension ExtAnalyze.

Standalone Usage:
-----------------

>>> pt = espresso.analysis.PressureTensor(system)
>>> print "pressure tensor of current configuration = ", pt.compute()

or 

>>> pt = espresso.analysis.PressureTensor(system)
>>> for k in range(100):
>>>     integrator.run(100)
>>>     pt.performMeasurement()
>>> print "average pressure tensor = ", pt.getAverageValue()

Usage in integrator with ExtAnalyze:
------------------------------------

>>> pt           = espresso.analysis.PressureTensor(system)
>>> extension_pt = espresso.integrator.ExtAnalyze(pt , interval=100)
>>> integrator.addExtension(extension_pt)
>>> integrator.run(10000)
>>> pt_ave = pt.getAverageValue()
>>> print "average Pressure Tensor = ", pt_ave[:6]
>>> print "          std deviation = ", pt_ave[6:]
>>> print "number of measurements  = ", pt.getNumberOfMeasurements()

The following methods are supported:

* performMeasurement()
    computes the pressure tensor and updates average and standard deviation
* reset()
    resets average and standard deviation to 0
* compute()
    computes the instant pressure tensor, return value: [xx, yy, zz, xy, xz, yz]
* getAverageValue()
    returns the average pressure tensor and the standard deviation,
    return value: [xx, yy, zz, xy, xz, yz, +-xx, +-yy, +-zz, +-xy, +-xz, +-yz]
* getNumberOfMeasurements()
    counts the number of measurements that have been computed (standalone or in integrator)
    does _not_ include measurements that have been done using "compute()"
"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.AnalysisBase import *
from _espresso import analysis_PressureTensor

class PressureTensorLocal(AnalysisBaseLocal, analysis_PressureTensor):
    'The (local) compute of pressure tensor.'
    def __init__(self, system):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_PressureTensor, system)

if pmi.isController:
    class PressureTensor(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.analysis.PressureTensorLocal',
            )
