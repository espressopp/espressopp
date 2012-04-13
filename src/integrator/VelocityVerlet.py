from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.MDIntegrator import *
from _espresso import integrator_VelocityVerlet 

class VelocityVerletLocal(MDIntegratorLocal, integrator_VelocityVerlet):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_VelocityVerlet, system)

if pmi.isController :
    class VelocityVerlet(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.VelocityVerletLocal',
            pmiproperty = [ 'langevin', 'berendsenBarostat', 'berendsenThermostat', 'langevinBarostat', 'isokinetic', 'stochasticVelocityRescaling'],
            pmicall = [ 'getTimers', 'resetTimers' ]
            )
