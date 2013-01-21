from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *
from _espresso import integrator_DPDThermostat 

class DPDThermostatLocal(ExtensionLocal, integrator_DPDThermostat):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_DPDThermostat, system, vl)

    #def enableAdress(self):
    #    if pmi.workerIsActive():
    #        self.cxxclass.enableAdress(self);

if pmi.isController :
    class DPDThermostat(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.DPDThermostatLocal',
            pmiproperty = [ 'gamma', 'tgamma', 'temperature' ]
            )
