"""
******************************************
**espresso.integrator.LangevinThermostat**
******************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *
from _espresso import integrator_LangevinThermostat 

class LangevinThermostatLocal(ExtensionLocal, integrator_LangevinThermostat):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LangevinThermostat, system)

    #def enableAdress(self):
    #    if pmi.workerIsActive():
    #        self.cxxclass.enableAdress(self);

if pmi.isController :
    class LangevinThermostat(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LangevinThermostatLocal',
            pmiproperty = [ 'gamma', 'temperature', 'adress' ]
            )
