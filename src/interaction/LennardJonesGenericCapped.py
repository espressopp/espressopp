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
********************************************
**espressopp.interaction.LennardJonesGenericCapped**
********************************************

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJonesGenericCapped, \
                      interaction_VerletListLennardJonesGenericCapped, \
                      interaction_VerletListAdressLennardJonesGenericCapped, \
                      interaction_VerletListHadressLennardJonesGenericCapped, \
                      interaction_CellListLennardJonesGenericCapped, \
                      interaction_FixedPairListLennardJonesGenericCapped

class LennardJonesGenericCappedLocal(PotentialLocal, interaction_LennardJonesGenericCapped):
    def __init__(self, epsilon=1.0, sigma=1.0, a=12, b=6,
                 cutoff=infinity, caprad=0.0, shift="auto"):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                print "here without shift"
                cxxinit(self, interaction_LennardJonesGenericCapped, 
                        epsilon, sigma, a, b, cutoff, caprad)
            else:
                print "here with shift"
                cxxinit(self, interaction_LennardJonesGenericCapped, 
                        epsilon, sigma, a, b, cutoff, caprad, shift)

class VerletListLennardJonesGenericCappedLocal(InteractionLocal, interaction_VerletListLennardJonesGenericCapped):
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJonesGenericCapped, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

    def getVerletListLocal(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getVerletList(self)

class VerletListAdressLennardJonesGenericCappedLocal(InteractionLocal, interaction_VerletListAdressLennardJonesGenericCapped):
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJonesGenericCapped, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
class VerletListHadressLennardJonesGenericCappedLocal(InteractionLocal, interaction_VerletListHadressLennardJonesGenericCapped):
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLennardJonesGenericCapped, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)
            
class CellListLennardJonesGenericCappedLocal(InteractionLocal, interaction_CellListLennardJonesGenericCapped):
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesGenericCapped, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesGenericCappedLocal(InteractionLocal, interaction_FixedPairListLennardJonesGenericCapped):
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJonesGenericCapped, system, vl, potential)
        
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
    class LennardJonesGenericCapped(Potential):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJonesGenericCappedLocal',
            pmiproperty = ['epsilon', 'sigma', 'a', 'b', 'cutoff', 'caprad']
            )

    class VerletListLennardJonesGenericCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListLennardJonesGenericCappedLocal',
            pmicall = ['setPotential', 'getPotential', 'getVerletList']
            )

    class VerletListAdressLennardJonesGenericCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLennardJonesGenericCappedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListAdressLennardJonesGenericCapped2(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListAdressLennardJonesGenericCapped2Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressLennardJonesGenericCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressLennardJonesGenericCappedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )
            
    class VerletListHadressLennardJonesGenericCapped2(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListHadressLennardJonesGenericCapped2Local',
            pmicall = ['setPotentialAT', 'setPotentialCG']
            )

    class CellListLennardJonesGenericCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListLennardJonesGenericCappedLocal',
            pmicall = ['setPotential']
            )
        
    class FixedPairListLennardJonesGenericCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListLennardJonesGenericCappedLocal',
            pmicall = ['getPotential', 'setPotential', 'setFixedPairList','getFixedPairList' ]
        )
