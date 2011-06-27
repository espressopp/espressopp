from espresso import MultiSystem
from espresso import pmi
import random
from math import exp

class ParallelTempering(object):
    
    def __init__(self, NumberOfSystems = 4, coupleEveryNsteps = 100):
        if (NumberOfSystems < 2):
            print "ERROR: Parallel Tempering needs at least 2 systems to be able to work"
        else:
            random.seed(12345)
            print "WARNING: random.seed set to 12345 ... still need to workout some more randomness"
            self._multisystem    = MultiSystem()
            self._nsystems       = NumberOfSystems
            self._couplesteps    = coupleEveryNsteps
            self._cpugroup       = []
            self._comm           = []
            self._ncpustotal     = pmi.size
            self._ncpuspersystem = self._ncpustotal / self._nsystems
            self._oddeven        = 0
            if (self._ncpuspersystem * self._nsystems != self._ncpustotal):
                print "ERROR: Number of Parallel Tempering systems times CPUs per system does not match total number of CPUs"
            else:
                for i in range(self._nsystems):
                    self._cpugroup.append( range(i * self._ncpuspersystem, (i+1) * self._ncpuspersystem) )
                    self._comm.append(pmi.Communicator(self._cpugroup[i]))
            
    def startDefiningSystem(self,n):
        if not (n in range(0,self._nsystems)):
            print "ERROR: System number must be between 0 and ",self._nsystems
        else:
            pmi.activate(self._comm[n])
            self._multisystem.beginSystemDefinition()

    def endDefiningSystem(self,n):
        if not (n in range(0,self._nsystems)):
            print "ERROR: System number must be between 0 and ",self._nsystems
        else:
            pmi.deactivate(self._comm[n])
            
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
        
    def setAnalysisNPart(self,analysisNPart):
        self._multisystem.setAnalysisNPart(analysisNPart)
        
    def run(self, nsteps):
        totalsteps = 0
        while totalsteps<nsteps:
          if self._oddeven == 0:
              self._oddeven = 1
          else:
              self._oddeven = 0;      
          self._multisystem.runIntegrator(self._couplesteps)
          energies     = self._multisystem.runAnalysisPotential()
          temperatures = self._multisystem.runAnalysisTemperature()
          for i in range(len(energies)/2):
              m = 2 * i + self._oddeven
              n = m + 1
              if n<len(energies):
                  metro = random.random()
                  t1    = temperatures[m]
                  t2    = temperatures[n]
                  e1    = temperatures[m]
                  e2    = temperatures[n]
                  delta = exp(-(e1-e2)*(1/t2 -1/t1))
                  if delta >= metro:
                      temperatures[n] = t1
                      temperatures[m] = t2
# TODO       self._multisystem.setTemperaturesScaled(temperatures)
          totalsteps += self._couplesteps
          
