from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *
from _espresso import integrator_LangevinThermostat1D 

class LangevinThermostat1DLocal(ExtensionLocal, integrator_LangevinThermostat1D):
    'The (local) Langevin Thermostat (1D).'
    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LangevinThermostat1D, system)

    #def enableAdress(self):
    #    if pmi.workerIsActive():
    #        self.cxxclass.enableAdress(self);

if pmi.isController :
    class LangevinThermostat1D(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LangevinThermostat1DLocal',
            pmiproperty = [ 'gamma', 'temperature', 'adress', 'direction' ]
            )
