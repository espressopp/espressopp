from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import integrator_Berendsen

class BerendsenLocal(integrator_Berendsen):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_Berendsen, system)

if pmi.isController :
    class Berendsen(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.BerendsenLocal',
            pmiproperty = [ 'tau', 'pressure' ]
            )
