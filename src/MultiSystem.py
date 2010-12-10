from espresso.esutil import cxxinit
from espresso import pmi
import MPI

class MultiSystemLocal(object):
    """Local MultiSystem to simulate and analyze several systems in parallel."""
    
    def __init__(self):
        pass

    def beginSystemDefinition(self):
        if pmi._PMIComm and pmi._PMIComm.isActive():
            if pmi._PMIComm.getMPIsubcomm() != MPI.COMM_NULL:
                self.groupRank = pmi._PMIComm.getMPIsubcomm().rank
            else:
                print "_MPIsubcomm is MPI.COMM_NULL"
        else:
            self.groupRank = pmi._MPIcomm.rank

    def setIntegrator(self, integrator):
        self.integrator = integrator

    def setAnalysisTemperature(self, temperature):
        self.analysisTemperature = temperature

    def setAnalysisPotential(self, potential):
        self.analysisPotential = potential

    def runIntegrator(self, niter):
        self.integrator.cxxclass.run(self.integrator, niter)

    def runAnalysisTemperature(self):
        if self.groupRank == 0:
            return self.analysisTemperature.cxxclass.compute(self.analysisTemperature)
        else :
            self.analysisTemperature.cxxclass.compute(self.analysisTemperature)

    def runAnalysisPotential(self):
        if self.groupRank == 0:
            return self.analysisPotential.cxxclass.computeEnergy(self.analysisPotential)
        else :
            self.analysisPotential.cxxclass.computeEnergy(self.analysisPotential)

if pmi.isController :
    class MultiSystem(object):
        """MultiSystemIntegrator to simulate and analyze several systems in parallel."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.MultiSystemLocal',
            pmicall = [ 'setIntegrator', 'runIntegrator', 'setAnalysisTemperature', 'beginSystemDefinition',
                        'setAnalysisPotential'],
            pmiinvoke = [ 'runAnalysisTemperature', 'runAnalysisPotential' ]
            )
