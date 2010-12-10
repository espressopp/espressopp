from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import integrator_Langevin 

class LangevinLocal(integrator_Langevin):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_Langevin, system)

if pmi.isController :
    class Langevin(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LangevinLocal',
            pmiproperty = [ 'gamma', 'temperature' ]
            )
