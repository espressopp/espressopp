from espresso import pmi

from espresso.integrator.MDIntegrator import *
from _espresso import integrator_VelocityVerlet

class VelocityVerletLocal(MDIntegratorLocal):
    def __init__(self, set, 
                 posProperty, velProperty, forceProperty,
                 timestep):
        if not hasattr(self, 'cxxobject'):
            self.cxxobject = \
                integrator_VelocityVerlet(set.cxxobject, 
                                          posProperty.cxxobject,
                                          velProperty.cxxobject, 
                                          forceProperty.cxxobject)

if pmi.IS_CONTROLLER:
    class VelocityVerlet(MDIntegrator):
        'The Velocity-Verlet integrator.'
        def __init__(self, set, 
                     posProperty, velProperty, forceProperty, 
                     timestep = 1):
            if not hasattr(self, 'pmiobject'):
                self.pmiobject = pmi.create("espresso.integrator.VelocityVerletLocal",
                                            set.pmiobject,
                                            posProperty.pmiobject,
                                            velProperty.pmiobject,
                                            forceProperty.pmiobject,
                                            timestep)
    
