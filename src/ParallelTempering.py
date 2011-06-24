from espresso import MultiSystem
from espresso import pmi
import random

class ParallelTempering(object):
    
    def __init__(self, NumberOfSystems = nsystems, coupleEveryNsteps=couplesteps):
        if (nsystems < 2):
            print "ERROR: Parallel Tempering needs at least 2 systems to be able to work"
        else:
            self._multisystem    = espresso.MultiSystem()
            self._nsystems       = nsystems
            self._couplesteps    = couplesteps
            self._cpugroup       = []
            self._comm           = []
            self._ncpustotal     = espresso.pmi.size
            self._ncpuspersystem = self._ncpustotal / self._nsystems
            if (self._ncpuspersystem * self._nsystem != self._ncpustotal):
                print "ERROR: Number of Parallel Tempering systems times CPUs per system does not match total number of CPUs"
            else:
                for i in range(self._nsystems):
                    self._cpugroup.append( range(i * self._ncpuspersystem, (i+1) * self._ncpuspersystem) )
                    self._comm.append(espresso.pmi.Communicator(self._cpugroup[i]))
            
    def startDefiningSystem(self,n):
        if not (n in range(0,self._nsystems)):
            print "ERROR: System number must be between 0 and ",self._nsystems
        else:
            espresso.pmi.activate(self._comm[n])
            self._multisystem.beginSystemDefinition()

    def endDefiningSystem(self,n):
        if not (n in range(0,self._nsystems)):
            print "ERROR: System number must be between 0 and ",self._nsystems
        else:
            espresso.pmi.deactivate(comm[n])
            
    def getNumberOfSystems(self):
        return self._nsystems
    
    def getNumberOfCPUsPerSystem(self):
        return self._ncpuspersystem
    
    def setIntegrator(self, integrator):
        self._multisystem.setIntegrator(integrator)
        
    def setAnalysisE(self, analysisE):
        self._multisystem.setAnalysisPotential(analysisE)
                
    def setAnalysisT(self, analysisT):
        self._multisystem.setAnalysisTemperature(analysisT)
        
    def run(selfself, nsteps):
        totalsteps = 0
        while totalsteps<nsteps:
          _multisystem.run(couplesteps)
          totalsteps += couplesteps
          
