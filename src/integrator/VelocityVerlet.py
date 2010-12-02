from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.MDIntegrator import *
from _espresso import integrator_VelocityVerlet 

class VelocityVerletLocal(MDIntegratorLocal, integrator_VelocityVerlet):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system):
        cxxinit(self, integrator_VelocityVerlet, system)

if pmi.isController :
    class VelocityVerlet(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.VelocityVerletLocal',
            pmiproperty = [ 'langevin' ],
            pmicall = [ 'getTimers' ]
            )
