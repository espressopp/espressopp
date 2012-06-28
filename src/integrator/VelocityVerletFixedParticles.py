from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.VelocityVerlet import *
from _espresso import integrator_VelocityVerletFixedParticles

class VelocityVerletFixedParticlesLocal(VelocityVerletLocal, integrator_VelocityVerletFixedParticles):
    'The (local) Velocity Verlet Integrator for fixed particles.'
    def __init__(self, system, fixedParticles, fixMask):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_VelocityVerletFixedParticles, system, fixedParticles, fixMask)

if pmi.isController :
    class VelocityVerletFixedParticles(VelocityVerlet):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.integrator.VelocityVerletFixedParticlesLocal',
          pmiproperty = [ 'fixedParticles', 'fixMask' ]
        )
