from espresso import pmi
from espresso.esutil import *

from espresso.integrator.MDIntegrator import *
from _espresso import integrator_VelocityVerlet

class VelocityVerletLocal(MDIntegratorLocal, integrator_VelocityVerlet):
    def __init__(self, set, posProperty, velProperty, forceProperty, timestep=0.1):
        cxxinit(self, integrator_VelocityVerlet, 
                set, posProperty, velProperty, forceProperty, timestep)

if pmi.IS_CONTROLLER:
    class VelocityVerlet(MDIntegrator):
        'The Velocity-Verlet integrator.'
#        pmiproxydefs = MDIntegrator.pmiproxydefs
#        pmiproxydefs['subjectclass'] = 'espresso.integrator.VelocityVerletLocal'
        
        def __init__(self, set, 
                     posProperty, velProperty, forceProperty,
                     timestep):
            pmiinit(self, "espresso.integrator.VelocityVerletLocal",
                    set.pmiobject,
                    posProperty.pmiobject,
                    velProperty.pmiobject,
                    forceProperty.pmiobject, 
                    timestep)    

        def run(self, steps):
            pmi.call(self.pmiobject.run, steps)
