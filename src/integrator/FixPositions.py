from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import integrator_FixPositions 

class FixPositionsLocal(integrator_FixPositions):
    'The (local) Fix Positions part.'
    def __init__(self, system, particleGroup, fixMask):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_FixPositions, system, particleGroup, fixMask)

if pmi.isController :
    class FixPositions(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.FixPositionsLocal',
            pmiproperty = [ 'particleGroup', 'fixMask' ]
            )
