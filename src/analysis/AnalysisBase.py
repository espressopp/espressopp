"""
***************************************************************
**AnalysisBase** - abstract base class for analysis/measurement
***************************************************************

This abstract base class provides the interface and some basic
functionality for classes that do analysis or observable measurements
  
It provides the following methods:

* performMeasurement()
    computes the observable and updates average and standard deviation
* reset()
    resets average and standard deviation
* compute()
    computes the instant value of the observable, return value is a python list or a scalar
* getAverageValue()
    returns the average value for the observable and the standard deviation,
    return value is a python list
* getNumberOfMeasurements()
    counts the number of measurements that have been performed (standalone or in integrator)
    does _not_ include measurements that have been done using "compute()"
"""

from espresso import pmi

from espresso.ParticleAccess import *
from _espresso import analysis_AnalysisBase

class AnalysisBaseLocal(ParticleAccessLocal, analysis_AnalysisBase):
    """Abstract local base class for observables."""
    def performMeasurement(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.performMeasurement(self)

    def reset(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.reset(self)

    def compute(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            res = self.cxxclass.compute(self)
            if len(res) > 1:
                return res
            else:
                return res[0]

    def getAverageValue(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getAverageValue(self)

    def getNumberOfMeasurements(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNumberOfMeasurements(self)

if pmi.isController :
    class AnalysisBase(ParticleAccess):
        """Abstract base class for observable."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "performMeasurement", "reset", "compute", "getAverageValue", "getNumberOfMeasurements" ]
            )
