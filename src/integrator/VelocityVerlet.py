from espresso import pmi

from espresso.integrator.MDIntegrator import *
from _espresso import integrator_VelocityVerlet

class VelocityVerletLocal(MDIntegratorLocal, integrator_VelocityVerlet):
    def __init__(self, set, posProperty, velProperty, forceProperty):
        if not hasattr(self, 'cxxinit'):
            integrator_VelocityVerlet.__init__(self,
                set, posProperty, velProperty, forceProperty)
            self.cxxinit = True

if pmi.IS_CONTROLLER:
    class VelocityVerlet(MDIntegrator):
        'The Velocity-Verlet integrator.'
        pmiproxydefs = MDIntegrator.pmiproxydefs
        pmiproxydefs['subjectclass'] = 'espresso.integrator.VelocityVerletLocal'
        
        def __init__(self, set, 
                     posProperty, velProperty, forceProperty, 
                     timestep=0.1):
            if not hasattr(self, 'pmiobject'):
                self.pmiobject = \
                    pmi.create("espresso.integrator.VelocityVerletLocal",
                               set.pmiobject,
                               posProperty.pmiobject,
                               velProperty.pmiobject,
                               forceProperty.pmiobject)
                self.timestep = timestep
    
