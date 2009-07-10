from espresso import pmi
import abc

class MDIntegratorLocal(object):
    __metaclass__ = abc.ABCMeta
    def step(self):
        return self.cxxobject.step()

    def run(self, steps = 1):
        return self.cxxobject.run(steps)

    @property
    def timeStep(self): return self.cxxobject.getTimeStep()
    @timeStep.setter
    def timeStep(self, _timeStep): self.cxxobject.setTimeStep(_timeStep)

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.integrator')
    class MDIntegrator(object):
        __metaclass__ = abc.ABCMeta
        pmiproxydefs = {
            'subjectclass': 'espresso.integrator.MDIntegrator',
            'pmicall' : [ 'run', 'step' ],
            'pmiproperty' : [ 'timeStep' ]
            }


    
