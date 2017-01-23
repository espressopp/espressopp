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
*************************************
espressopp.interaction.HarmonicUnique
*************************************

.. math::

        U = K  (d - d_{cur})^2;






.. function:: espressopp.interaction.HarmonicUnique(K)

		:param K: (default: 1.0)
		:type K: real

.. function:: espressopp.interaction.FixedPairDistListHarmonicUnique(system, fpl, potential)

		:param system: 
		:param fpl: 
		:param potential: 
		:type system: 
		:type fpl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairDistListHarmonicUnique.getFixedPairList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedPairDistListHarmonicUnique.setFixedPairList(fixedpairlist)

		:param fixedpairlist: 
		:type fixedpairlist: 

.. function:: espressopp.interaction.FixedPairDistListHarmonicUnique.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.PotentialUniqueDist import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_HarmonicUnique, \
                      interaction_FixedPairDistListHarmonicUnique

class HarmonicUniqueLocal(PotentialUniqueDistLocal, interaction_HarmonicUnique):

    def __init__(self, K=1.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, interaction_HarmonicUnique, K)

class FixedPairDistListHarmonicUniqueLocal(InteractionLocal, interaction_FixedPairDistListHarmonicUnique):

    def __init__(self, system, fpl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairDistListHarmonicUnique, system, fpl, potential)

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
    class HarmonicUnique(PotentialUniqueDist):
        'The HarmonicUnique potential.'
        pmiproxydefs = dict(
          cls = 'espressopp.interaction.HarmonicUniqueLocal',
          pmiproperty = ['K']
        )

    class FixedPairDistListHarmonicUnique(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.FixedPairDistListHarmonicUniqueLocal',
          pmicall = ['setPotential','setFixedPairList','getFixedPairList']
        )
