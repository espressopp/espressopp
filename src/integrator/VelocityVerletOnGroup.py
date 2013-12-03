"""
*********************************************
**espresso.integrator.VelocityVerletOnGroup**
*********************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.MDIntegrator import *
from _espresso import integrator_VelocityVerletOnGroup

class VelocityVerletOnGroupLocal(MDIntegratorLocal, integrator_VelocityVerletOnGroup):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system, group):
        cxxinit(self, integrator_VelocityVerletOnGroup, system, group)

if pmi.isController :
    class VelocityVerletOnGroup(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.VelocityVerletOnGroupLocal',
            pmiproperty = [ 'langevin' ]
            )
