"""
*************************************
**ExtAnalyze** - Integrator Extension
*************************************

This class can be used to execute nearly all analysis objects
within the main integration loop which allows to automatically
accumulate time averages (with standard deviation error bars). 
  
Example Usage:
-----------------

>>> pt           = espresso.analysis.PressureTensor(system)
>>> extension_pt = espresso.integrator.ExtAnalyze(pt , interval=100)
>>> integrator.addExtension(extension_pt)
>>> integrator.run(10000)
>>>
>>> pt_ave = pt.getAverageValue()
>>> print "average Pressure Tensor = ", pt_ave[:6]
>>> print "          std deviation = ", pt_ave[6:]
>>> print "number of measurements  = ", pt.getNumberOfMeasurements()
"""

from espresso.esutil import cxxinit
from espresso import pmi
from espresso.integrator.Extension import *
from _espresso import integrator_ExtAnalyze 

class ExtAnalyzeLocal(ExtensionLocal, integrator_ExtAnalyze):
    'The (local) extension analyze.'
    def __init__(self, action_obj, interval=1):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
           cxxinit(self, integrator_ExtAnalyze, action_obj, interval)

if pmi.isController :
    class ExtAnalyze(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.ExtAnalyzeLocal',
        )
