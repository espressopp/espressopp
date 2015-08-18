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
***********************************
**espressopp.interaction.Harmonic**
***********************************

.. math::
	U = K (d - r_0)^2


.. function:: espressopp.interaction.Harmonic(K, r0, cutoff, shift)

		:param K: (default: 1.0)
		:param r0: (default: 0.0)
		:param cutoff: (default: infinity)
		:param shift: (default: 0.0)
		:type K: real
		:type r0: real
		:type cutoff: 
		:type shift: real

.. function:: espressopp.interaction.FixedPairListHarmonic(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListHarmonic.getFixedPairList()

		:rtype: FixedPairList (A python list of lists, *not* tuples)

.. function:: espressopp.interaction.FixedPairListHarmonic.setFixedPairList(fixedpairlist)

		:param fixedpairlist: 
		:type fixedpairlist: FixedPairList


.. function:: espressopp.interaction.FixedPairListHarmonic.setPotential(potential)

		:param potential: 
		:type potential: 

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_Harmonic, interaction_FixedPairListHarmonic

class HarmonicLocal(PotentialLocal, interaction_Harmonic):

    def __init__(self, K=1.0, r0=0.0, 
                 cutoff=infinity, shift=0.0):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift == "auto":
                cxxinit(self, interaction_Harmonic, K, r0, cutoff)
            else:
                cxxinit(self, interaction_Harmonic, K, r0, cutoff, shift)

class FixedPairListHarmonicLocal(InteractionLocal, interaction_FixedPairListHarmonic):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListHarmonic, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    
    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

if pmi.isController:
    class Harmonic(Potential):

        pmiproxydefs = dict(
            cls = 'espressopp.interaction.HarmonicLocal',
            pmiproperty = ['K', 'r0']
            )

    class FixedPairListHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListHarmonicLocal',
            pmicall = ['setPotential','getPotential','setFixedPairList','getFixedPairList']
            )
