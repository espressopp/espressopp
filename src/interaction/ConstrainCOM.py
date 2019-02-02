#  Copyright (C) 2017, 2019
#      Max Planck Institute for Polymer Research
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
espressopp.interaction.ConstrainCOM
***********************************

This class is for calculating forces acting on constrained center of mass of subchains [Zhang_2014]_.

Subchains are defined as a tuple list.

.. math:: U = k_{com} \left(\vec{r_{com}} - \vec{R_{com}}\right)^2,

where :math:`\vec{r_{com}}` stands for the center of mass of subchain and :math:`\vec{R_{com}}`
stands for the desired center of mass of subchain.

This class implies 2 conditions on a tuple list defining subchains:

1. The length of all tuples must be the same.
2. int(key particle id / The length of a tuple) must not be redundantly, where key particle id
   is the smallest particle id in a tuple.


.. [Zhang_2014] G. Zhang, L. A. Moreira, T. Stuehn, K. Ch. Daoulas, and K. Kremer,
         Equilibration of high molecular weight polymer melts: A hierarchical strategy, Macro Lett., 2014, 3, 198

.. function:: espressopp.interaction.ConstrainCOM(k_com)

		:param k_com: (default: 100.)
		:type k_com: real

.. function:: espressopp.interaction.FixedLocalTupleListConstrainCOM(system, tuplelist, potential)

		:param system: 
		:param tuplelist: 
		:param potential: 
		:type system: 
		:type tuplelist: 
		:type potential: 

.. function:: espressopp.interaction.FixedLocalTupleListConstrainCOM.getPotential()

		:rtype: 

.. function:: espressopp.interaction.FixedLocalTupleListConstrainCOM.setCom(particlelist)

		:param particlelist: 
		:type particlelist:


"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_ConstrainCOM, interaction_FixedLocalTupleListConstrainCOM

class ConstrainCOMLocal(PotentialLocal, interaction_ConstrainCOM):
    
    def __init__(self, k_com=100.):
        """Initialize the local ConstrainCOM."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_ConstrainCOM, k_com)
        
class FixedLocalTupleListConstrainCOMLocal(InteractionLocal, interaction_FixedLocalTupleListConstrainCOM):
    
    def __init__(self, system, fixedtuplelist, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedLocalTupleListConstrainCOM, system, fixedtuplelist, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

    def setCom(self, particleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            id = 0
            for particle in particleList:
                pos = particle[0]
                self.cxxclass.setCom(self, id, pos)
                id = id + 1

if pmi.isController:
    class ConstrainCOM(Potential):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.ConstrainCOMLocal',
            pmiproperty = ['k_com'],
            )

    class FixedLocalTupleListConstrainCOM(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedLocalTupleListConstrainCOMLocal',
            pmicall = ['getPotential', 'setCom']
            )

