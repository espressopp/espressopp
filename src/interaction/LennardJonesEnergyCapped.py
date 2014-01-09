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


"""
*************************************************
**espresso.interaction.LennardJonesEnergyCapped**
*************************************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_LennardJonesEnergyCapped, \
                      interaction_VerletListLennardJonesEnergyCapped, \
                      interaction_VerletListAdressLennardJonesEnergyCapped, \
                      interaction_VerletListHadressLennardJonesEnergyCapped, \
                      interaction_CellListLennardJonesEnergyCapped, \
                      interaction_FixedPairListLennardJonesEnergyCapped

class LennardJonesEnergyCappedLocal(PotentialLocal, interaction_LennardJonesEnergyCapped):
    'The (local) Lennard-Jones potential with energy capping.'
    def __init__(self, epsilon=1.0, sigma=1.0, 
                 cutoff=infinity, caprad=0.0 ,shift="auto"):
        """Initialize the local Lennard Jones object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_LennardJonesEnergyCapped, 
                        epsilon, sigma, cutoff, caprad)
            else:
                cxxinit(self, interaction_LennardJonesEnergyCapped, 
                        epsilon, sigma, cutoff, caprad, shift)

class VerletListLennardJonesEnergyCappedLocal(InteractionLocal, interaction_VerletListLennardJonesEnergyCapped):
    'The (local) Lennard Jones interaction using Verlet lists.'
    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJonesEnergyCapped, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class VerletListAdressLennardJonesEnergyCappedLocal(InteractionLocal, interaction_VerletListAdressLennardJonesEnergyCapped):
    'The (local) Lennard Jones interaction using Verlet lists.'
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListAdressLennardJonesEnergyCapped, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

    def getPotentialAT(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotentialAT(self, type1, type2)

    def getPotentialCG(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotentialCG(self, type1, type2)
        
class VerletListHadressLennardJonesEnergyCappedLocal(InteractionLocal, interaction_VerletListHadressLennardJonesEnergyCapped):
    'The (local) Lennard Jones interaction using Verlet lists.'
    def __init__(self, vl, fixedtupleList):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListHadressLennardJonesEnergyCapped, vl, fixedtupleList)

    def setPotentialAT(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialAT(self, type1, type2, potential)

    def setPotentialCG(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotentialCG(self, type1, type2, potential)

    def getPotentialAT(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotentialAT(self, type1, type2)

    def getPotentialCG(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotentialCG(self, type1, type2)

class CellListLennardJonesEnergyCappedLocal(InteractionLocal, interaction_CellListLennardJonesEnergyCapped):
    'The (local) Lennard Jones interaction using cell lists.'
    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesEnergyCapped, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class FixedPairListLennardJonesEnergyCappedLocal(InteractionLocal, interaction_FixedPairListLennardJonesEnergyCapped):
    'The (local) Lennard-Jones interaction using FixedPair lists.'
    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJonesEnergyCapped, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getPotential(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self)

if pmi.isController:
    class LennardJonesEnergyCapped(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
            cls = 'espresso.interaction.LennardJonesEnergyCappedLocal',
            pmiproperty = ['epsilon', 'sigma', 'cutoff', 'caprad']
            )

    class VerletListLennardJonesEnergyCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListLennardJonesEnergyCappedLocal',
            pmicall = ['setPotential', 'getPotential']
            )

    class VerletListAdressLennardJonesEnergyCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListAdressLennardJonesEnergyCappedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG', 'getPotentialAT', 'getPotentialCG']
            )
            
    class VerletListHadressLennardJonesEnergyCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.VerletListHadressLennardJonesEnergyCappedLocal',
            pmicall = ['setPotentialAT', 'setPotentialCG', 'getPotentialAT', 'getPotentialCG']
            )

    class CellListLennardJonesEnergyCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.CellListLennardJonesEnergyCappedLocal',
            pmicall = ['setPotential', 'getPotential']
            )
        
    class FixedPairListLennardJonesEnergyCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.interaction.FixedPairListLennardJonesEnergyCappedLocal',
            pmicall = ['setPotential']
            )
