from espresso import esutil, pmi

from _espresso import integrator_VelocityVerlet as _VelocityVerlet

class VelocityVerletLocal(_VelocityVerlet):
    def __init__(self, particles, 
                 positionProperty, velocityProperty, forceProperty, 
                 timeStep):
        # TODO: defaults:
        # - if "velocity", "position", "force" exist in global table, use them
        # - if not, create them
        _VelocityVerlet.__init__(self, particles, 
                                 positionProperty, 
                                 velocityProperty, 
                                 forceProperty)
        self.setTimeStep(timeStep)

    def run(self, steps = 1):
        return _VelocityVerlet.run(self, steps)

    @property
    def timeStep(self): return self.getTimeStep
    @timeStep.setter
    def timeStep(self, _timeStep): self.setTimeStep(_timeStep)

    if pmi.IS_CONTROLLER:
        pmi.exec_('import espresso.integrator.VelocityVerlet')

        class VelocityVerlet (object):
            'The Velocity-Verlet integrator.'
            __metaclass__ = pmi.Proxy
            pmiproxydefs = {
                'subjectclass': 'espresso.integrator.VelocityVerletLocal',
                'pmicall' : [ 'run' ],
                'pmiproperty' : [ 'timeStep' ]
                }

    
