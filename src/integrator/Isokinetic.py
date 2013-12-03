"""
**********************************
**espresso.integrator.Isokinetic**
**********************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *
from _espresso import integrator_Isokinetic 

class IsokineticLocal(ExtensionLocal, integrator_Isokinetic):
    'The (local) Isokinetic Thermostat.'
    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_Isokinetic, system)

if pmi.isController :
    class Isokinetic(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.integrator.IsokineticLocal',
          pmiproperty = [ 'temperature', 'coupling' ]
        )
