from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.MDIntegrator import *
from _espresso import integrator_VelocityVerletAdress 

class VelocityVerletAdressLocal(MDIntegratorLocal, integrator_VelocityVerletAdress):
    'The (local) Velocity Verlet Adress Integrator.'
    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_VelocityVerletAdress, system)

if pmi.isController :
    class VelocityVerletAdress(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.VelocityVerletAdressLocal',
            pmiproperty = [ 'langevin' ],
            pmicall = [ 'getTimers' ]
            )
