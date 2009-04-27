from espresso import esutil
from espresso.esutil import choose
from espresso import pmi

# wrap VelocityVerlet
pmi.exec_('from _espresso import integrator_VelocityVerlet')
class VelocityVerlet (object):
    'The velocity-Verlet integrator.'

    def __init__(self, timeStep) :
        """Initialize the (parallel) Velocity Verlet object."""
        object.__init__(self)
        # create the pmi object
        self.worker = pmi.create('integrator_VelocityVerlet', timeStep)
        # set the defaults
        #self.set(epsilon, sigma, cutoff)

    

