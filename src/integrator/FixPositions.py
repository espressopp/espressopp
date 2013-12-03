"""
************************************
**espresso.integrator.FixPositions**
************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi
from espresso.integrator.Extension import *
from _espresso import integrator_FixPositions 

class FixPositionsLocal(ExtensionLocal, integrator_FixPositions):
    'The (local) Fix Positions part.'
    def __init__(self, system, particleGroup, fixMask):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_FixPositions, system, particleGroup, fixMask)

if pmi.isController :
    class FixPositions(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.FixPositionsLocal',
            pmicall = ['setFixMask', 'getFixMask'],
            pmiproperty = [ 'particleGroup' ]
            )
