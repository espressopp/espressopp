# Copyright (C) 2012-2016 Max Planck Institute for Polymer Research
# Copyright (C) 2008-2011 Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
espressopp.interaction.LJcos
****************************

if :math:`r^2 \leq border_{pot}`, then:

.. math::

	U = 4(\frac{1}{r^{12}} - \frac{1}{r^6}) + 1 - \phi

else:
	
.. math::     

	U = \frac{1}{2}\phi (cos(\alpha r^2+\beta) - 1)







.. function:: espressopp.interaction.LJcos(phi)

		:param phi: (default: 1.0)
		:type phi: real

.. function:: espressopp.interaction.VerletListLJcos(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListLJcos.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListLJcos.getVerletList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.VerletListLJcos.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressLJcos(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListAdressLJcos.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressLJcos.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressLJcos(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListHadressLJcos.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressLJcos.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListLJcos(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListLJcos.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLJcos(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLJcos.getFixedPairList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedPairListLJcos.setFixedPairList(fixedpairlist)

		:param fixedpairlist: 
		:type fixedpairlist: 

.. function:: espressopp.interaction.FixedPairListLJcos.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LJcos, \
                      interaction_VerletListLJcos, \
                      interaction_VerletListAdressLJcos, \
                      interaction_VerletListHadressLJcos, \
                      interaction_CellListLJcos, \
                      interaction_FixedPairListLJcos

class LJcosLocal(PotentialLocal, interaction_LJcos):

    def __init__(self, phi=1.0):

      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, interaction_LJcos, phi)

class VerletListLJcosLocal(InteractionLocal, interaction_VerletListLJcos):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLJcos, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressLJcosLocal(InteractionLocal, interaction_VerletListAdressLJcos):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLJcos, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
class VerletListHadressLJcosLocal(InteractionLocal, interaction_VerletListHadressLJcos):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLJcos, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class CellListLJcosLocal(InteractionLocal, interaction_CellListLJcos):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLJcos, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLJcosLocal(InteractionLocal, interaction_FixedPairListLJcos):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLJcos, system, vl, potential)
        
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
    class LJcos(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LJcosLocal',
            pmiproperty = ['phi','sigma']
        )

    class VerletListLJcos(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListLJcosLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressLJcos(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLJcosLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressLJcos(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressLJcosLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListLJcos(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListLJcosLocal',
            pmicall = ['setPotential']
            )
        
    class FixedPairListLJcos(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListLJcosLocal',
            pmicall = ['setPotential', 'setFixedPairList','getFixedPairList' ]
        )
