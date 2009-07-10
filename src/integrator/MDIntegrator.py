from espresso import pmi
import abc

class MDIntegratorLocal(object):
    __metaclass__ = abc.ABCMeta
    def step(self):
        return self.cxxobject.step()
    def run(self, steps = 1):
        return self.cxxobject.run(steps)
    @property
    def timestep(self): return self.cxxobject.getTimeStep()
    @timestep.setter
    def timestep(self, _timeStep): self.cxxobject.setTimeStep(_timeStep)

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.integrator')
    class MDIntegrator(object):
        __metaclass__ = abc.ABCMeta
        def run(self, steps = 1):
            return pmi.call(self.pmiobject.run, steps)

        def step(self):
            return pmi.call(self.pmiobject.step)

        @property
        def timestep(self): 
            return self.pmiobject.timestep

        @timestep.setter
        def timestep(self, _timestep): 
            pmi.call("espresso.integrator.MDIntegratorLocal.timestep.fset", self.pmiobject, _timestep)


    
