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
**********************
espressopp.MultiSystem
**********************

.. function:: espressopp.MultiSystem()
.. function:: espressopp.MultiSystem.beginSystemDefinition()

		:rtype: 
		
.. function:: espressopp.MultiSystem.runAnalysisNPart()

		:rtype: 
		
.. function:: espressopp.MultiSystem.runAnalysisPotential()

		:rtype: 
		
.. function:: espressopp.MultiSystem.runAnalysisTemperature()

		:rtype: 
		
.. function:: espressopp.MultiSystem.runIntegrator(niter)

		:param niter: 
		:type niter: 
		:rtype: 
		
.. function:: espressopp.MultiSystem.setAnalysisNPart(npart)

		:param npart: 
		:type npart: 
		
.. function:: espressopp.MultiSystem.setAnalysisPotential(potential)

		:param potential: 
		:type potential: 
		
.. function:: espressopp.MultiSystem.setAnalysisTemperature(temperature)

		:param temperature: 
		:type temperature: 
		
.. function:: espressopp.MultiSystem.setIntegrator(integrator)

		:param integrator: 
		:type integrator: 
		
"""

from espressopp.esutil import cxxinit
from espressopp import pmi
import mpi4py.MPI as MPI

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

    def setAnalysisNPart(self, npart):
        self.analysisNPart = npart

    def setDumpConfXYZ(self, dumpconf):
        self.dumpConfXYZ = dumpconf

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

    def runAnalysisNPart(self):
        if self.groupRank == 0:
            return int(self.analysisNPart.cxxclass.compute(self.analysisNPart))
        else :
            self.analysisNPart.cxxclass.compute(self.analysisNPart)

    def runDumpConfXYZ(self):
        if self.groupRank == 0:
            return self.dumpConfXYZ.cxxclass.dump(self.dumpConfXYZ)
        else :
            self.dumpConfXYZ.cxxclass.dump(self.dumpConfXYZ)

if pmi.isController :
    class MultiSystem(object):
        """MultiSystemIntegrator to simulate and analyze several systems in parallel."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.MultiSystemLocal',
            pmicall = [ 'setIntegrator', 'runIntegrator', 'setAnalysisTemperature', 'beginSystemDefinition',
                        'setAnalysisPotential','setAnalysisNPart', 'setDumpConfXYZ'],
            pmiinvoke = [ 'runAnalysisTemperature', 'runAnalysisPotential','runAnalysisNPart','runDumpConfXYZ' ]
            )
