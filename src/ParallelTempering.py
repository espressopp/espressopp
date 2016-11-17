#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


r"""
****************************
espressopp.ParallelTempering
****************************

.. function:: espressopp.ParallelTempering(NumberOfSystems, RNG)

		:param NumberOfSystems: (default: 4)
		:param RNG: (default: None)
		:type NumberOfSystems: int
		:type RNG: 
		
.. function:: espressopp.ParallelTempering.endDefiningSystem(n)

		:param n: 
		:type n: 
		:rtype: 
		
.. function:: espressopp.ParallelTempering.exchange()

		:rtype: 
		
.. function:: espressopp.ParallelTempering.getNumberOfCPUsPerSystem()

		:rtype: 
		
.. function:: espressopp.ParallelTempering.getNumberOfSystems()

		:rtype: 
		
.. function:: espressopp.ParallelTempering.run(nsteps)

		:param nsteps: 
		:type nsteps: 
		:rtype: 
		
.. function:: espressopp.ParallelTempering.setAnalysisE(analysisE)

		:param analysisE: 
		:type analysisE: 
		
.. function:: espressopp.ParallelTempering.setAnalysisNPart(analysisNPart)

		:param analysisNPart: 
		:type analysisNPart: 
		
.. function:: espressopp.ParallelTempering.setAnalysisT(analysisT)

		:param analysisT: 
		:type analysisT: 
		
.. function:: espressopp.ParallelTempering.setIntegrator(integrator, thermostat)

		:param integrator: 
		:param thermostat: 
		:type integrator: 
		:type thermostat: 
		
.. function:: espressopp.ParallelTempering.startDefiningSystem(n)

		:param n: 
		:type n: 
		:rtype: 
"""

from espressopp import MultiSystem
from espressopp import pmi
import random
from math import exp

class ParallelTempering(object):
    
    def __init__(self, NumberOfSystems = 4, RNG = None):
        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        print "+ Due to a recent change in the design of ESPResSo++ the Parallel Tempering     +"
        print "+ Class is available in the current version.                                    +"
        print "+ Multisystem simulations are still possible but have to be setup manually.     +"
        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        #return 0 
        if (RNG == None):
            print "ERROR: ParallelTempering needs a random number generator"
        else: 
          if (NumberOfSystems < 2):
              print "ERROR: Parallel Tempering needs at least 2 systems to be able to work"
          else:
              self._RNG            = RNG
              self._multisystem    = MultiSystem()
              self._nsystems       = NumberOfSystems
              self._cpugroup       = []
              self._comm           = []
              self._thermostat     = []
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
    
    def setIntegrator(self, integrator, thermostat):
        self._multisystem.setIntegrator(integrator)
        self._thermostat.append(thermostat)
        
    def setAnalysisE(self, analysisE):
        self._multisystem.setAnalysisPotential(analysisE)
                
    def setAnalysisT(self, analysisT):
        self._multisystem.setAnalysisTemperature(analysisT)
        
    def setAnalysisNPart(self,analysisNPart):
        self._multisystem.setAnalysisNPart(analysisNPart)

    def setDumpConfXYZ(self, dumpconf):
        self._multisystem.setDumpConfXYZ(dumpconf)

    def runDumpConfXYZ(self):
        self._multisystem.runDumpConfXYZ()

    def run(self, nsteps):
        self._multisystem.runIntegrator(nsteps)

    def exchange(self):
        if self._oddeven == 0:
            self._oddeven = 1
        else:
            self._oddeven = 0;      
        energies     = self._multisystem.runAnalysisPotential()
        temperatures = self._multisystem.runAnalysisTemperature()
        nparts       = self._multisystem.runAnalysisNPart()
        print "energies     = ", energies
        print "temperatures = ", temperatures
        for i in range(len(energies)/2):
            m = 2 * i + self._oddeven
            n = m + 1
            if n<len(energies):
                metro = self._RNG.random()
                t1    = temperatures[m]
                t2    = temperatures[n]
                e1    = energies[m]
                e2    = energies[n]
                n1    = nparts[m]
                n2    = nparts[n]
                delta = (1/t2-1/t1)*(e1/n1-e2/n2)
                if delta <= 0:
                  exyesno = 'yes'
                else:
                  if exp(-delta) > metro:
                    exyesno = 'yes'
                  else:
                    exyesno = 'no'
                print "systems %i and %i: dE=%10.5f random=%10.5f ==> exchange: %s" % (m, n, delta, metro, exyesno)
                if exyesno=='yes':
                    # exchange temperature of system[n] <--> system[m]
                    pmi.activate(self._comm[n])
                    self._multisystem.beginSystemDefinition()
                    self._thermostat[n].temperature = t1
                    pmi.deactivate(self._comm[n])
                    pmi.activate(self._comm[m])
                    self._multisystem.beginSystemDefinition()
                    self._thermostat[m].temperature = t2
                    pmi.deactivate(self._comm[m])
          
