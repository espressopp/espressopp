#  Copyright (C) 2017
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
**********************************
espressopp.interaction.ConstrainRG
**********************************

This class calculates forces acting on constrained radii of gyration of subchains [Zhang_2014]_.

Subchains are defined as a tuple list.

.. math:: U = k_{rg} \left(R_{g}^2 - {R_{g}^{ideal}}^2\right)^2

where :math:`R_{g}^{ideal}` stands for the desired radius of gyration of subchain.

This class set 2 conditions on a tuple list. defining subchains.

1. The length of all tuples must be the same.
2. int(key particle id / The length of a tuple) must not be redundantly, where key particle id is the smallest particle id in a tuple.

.. function:: espressopp.interaction.ConstrainRG(k_rg)

		:param k_rg: (default: 100.)
		:type k_rg: real

.. function:: espressopp.interaction.FixedLocalTupleListConstrainRG(system, tuplelist, potential)

		:param system: 
		:param tuplelist: 
		:param potential: 
		:type system: 
		:type tuplelist: 
		:type potential: 

.. function:: espressopp.interaction.FixedLocalTupleListConstrainRG.getPotential()

		:rtype: 

.. function:: espressopp.interaction.FixedLocalTupleListConstrainRG.setRG(particlelist)

		:param particlelist:
		:type particlelist: python::list


"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_ConstrainRG, interaction_FixedLocalTupleListConstrainRG

class ConstrainRGLocal(PotentialLocal, interaction_ConstrainRG):
    
    def __init__(self, k_rg=100.):
        """Initialize the local ConstrainRG."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_ConstrainRG, k_rg)
        
class FixedLocalTupleListConstrainRGLocal(InteractionLocal, interaction_FixedLocalTupleListConstrainRG):
    
    def __init__(self, system, fixedtuplelist, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedLocalTupleListConstrainRG, system, fixedtuplelist, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

    def setRG(self, particleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            id = 0
            for particle in particleList:
                rg = particle[1]**2
                self.cxxclass.setRG(self, id, rg)
                id = id + 1

if pmi.isController:
    class ConstrainRG(Potential):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.ConstrainRGLocal',
            pmiproperty = ['k_rg'],
            )

    class FixedLocalTupleListConstrainRG(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedLocalTupleListConstrainRGLocal',
            pmicall = ['getPotential', 'setRG']
            )

