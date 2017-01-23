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
******************************************
espressopp.interaction.LennardJonesGeneric
******************************************

This class provides methods to compute forces and energies of
a generic Lennard Jones potential with arbitrary integers a and b.

.. math::

	V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{a} -
	\left( \frac{\sigma}{r} \right)^{b} \right]






.. function:: espressopp.interaction.LennardJonesGeneric(epsilon, sigma, a, b, cutoff, shift)

		:param epsilon: (default: 1.0)
		:param sigma: (default: 1.0)
		:param a: (default: 12)
		:param b: (default: 6)
		:param cutoff: (default: infinity)
		:param shift: (default: "auto")
		:type epsilon: real
		:type sigma: real
		:type a: int
		:type b: int
		:type cutoff: 
		:type shift: 

.. function:: espressopp.interaction.VerletListLennardJonesGeneric(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListLennardJonesGeneric.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListLennardJonesGeneric.getVerletList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.VerletListLennardJonesGeneric.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressLennardJonesGeneric(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListAdressLennardJonesGeneric.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressLennardJonesGeneric.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressLennardJonesGeneric2(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListAdressLennardJonesGeneric2.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListAdressLennardJonesGeneric2.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressLennardJonesGeneric(vl, fixedtupleList)

		:param vl: 
		:param fixedtupleList: 
		:type vl: 
		:type fixedtupleList: 

.. function:: espressopp.interaction.VerletListHadressLennardJonesGeneric.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressLennardJonesGeneric.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressLennardJonesGeneric2(vl, fixedtupleList, KTI)

		:param vl: 
		:param fixedtupleList: 
		:param KTI: (default: False)
		:type vl: 
		:type fixedtupleList: 
		:type KTI: 

.. function:: espressopp.interaction.VerletListHadressLennardJonesGeneric2.setPotentialAT(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.VerletListHadressLennardJonesGeneric2.setPotentialCG(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListLennardJonesGeneric(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListLennardJonesGeneric.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLennardJonesGeneric(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLennardJonesGeneric.getFixedPairList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedPairListLennardJonesGeneric.getPotential()

		:rtype: 

.. function:: espressopp.interaction.FixedPairListLennardJonesGeneric.setFixedPairList(fixedpairlist)

		:param fixedpairlist: 
		:type fixedpairlist: 

.. function:: espressopp.interaction.FixedPairListLennardJonesGeneric.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJonesGeneric, \
                      interaction_VerletListLennardJonesGeneric, \
                      interaction_VerletListAdressLennardJonesGeneric, \
                      interaction_VerletListAdressLennardJonesGeneric2, \
                      interaction_VerletListHadressLennardJonesGeneric, \
                      interaction_VerletListHadressLennardJonesGeneric2, \
                      interaction_CellListLennardJonesGeneric, \
                      interaction_FixedPairListLennardJonesGeneric

class LennardJonesGenericLocal(PotentialLocal, interaction_LennardJonesGeneric):

    def __init__(self, epsilon=1.0, sigma=1.0, a=12, b=6,
                 cutoff=infinity, shift="auto"):
        """Initialize the local generic Lennard Jones object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                print "here without shift"
                cxxinit(self, interaction_LennardJonesGeneric, 
                        epsilon, sigma, a, b, cutoff)
            else:
                print "here with shift"
                cxxinit(self, interaction_LennardJonesGeneric, 
                        epsilon, sigma, a, b, cutoff, shift)

class VerletListLennardJonesGenericLocal(InteractionLocal, interaction_VerletListLennardJonesGeneric):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJonesGeneric, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressLennardJonesGenericLocal(InteractionLocal, interaction_VerletListAdressLennardJonesGeneric):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJonesGeneric, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
class VerletListAdressLennardJonesGeneric2Local(InteractionLocal, interaction_VerletListAdressLennardJonesGeneric2):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJonesGeneric2, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
                  
class VerletListHadressLennardJonesGenericLocal(InteractionLocal, interaction_VerletListHadressLennardJonesGeneric):

    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLennardJonesGeneric, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
class VerletListHadressLennardJonesGeneric2Local(InteractionLocal, interaction_VerletListHadressLennardJonesGeneric2):

    def __init__(self, vl, fixedtupleList, KTI = False):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLennardJonesGeneric2, vl, fixedtupleList, KTI)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

class CellListLennardJonesGenericLocal(InteractionLocal, interaction_CellListLennardJonesGeneric):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesGeneric, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesGenericLocal(InteractionLocal, interaction_FixedPairListLennardJonesGeneric):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJonesGeneric, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

    def setFixedPairList(self, fixedpairlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setFixedPairList(self, fixedpairlist)

    
    def getFixedPairList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedPairList(self)

if pmi.isController:
    class LennardJonesGeneric(Potential):
        'The generic Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJonesGenericLocal',
            pmiproperty = ['epsilon', 'sigma', 'a', 'b']
            )

    class VerletListLennardJonesGeneric(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListLennardJonesGenericLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressLennardJonesGeneric(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLennardJonesGenericLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListAdressLennardJonesGeneric2(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLennardJonesGeneric2Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressLennardJonesGeneric(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressLennardJonesGenericLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressLennardJonesGeneric2(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressLennardJonesGeneric2Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListLennardJonesGeneric(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListLennardJonesGenericLocal',
            pmicall = ['setPotential']
            )
        
    class FixedPairListLennardJonesGeneric(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListLennardJonesGenericLocal',
            pmicall = ['getPotential', 'setPotential', 'setFixedPairList','getFixedPairList' ]
        )
