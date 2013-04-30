from espresso.esutil import cxxinit
from espresso import pmi
from espresso.integrator.Extension import *
from _espresso import integrator_LatticeBoltzmann 

class LatticeBoltzmannLocal(ExtensionLocal, integrator_LatticeBoltzmann):
    'The (local) Lattice Boltzmann part.'
    def __init__(self, system, particleGroup, fixMask):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LatticeBoltzmann, system, particleGroup, fixMask)

if pmi.isController :
    class LatticeBoltzmann(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LatticeBoltzmannLocal',
            pmicall = ['setFixMask', 'getFixMask'],
            pmiproperty = [ 'particleGroup' ]
            )
