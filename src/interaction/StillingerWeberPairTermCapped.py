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
******************************************************
**espressopp.interaction.StillingerWeberPairTermCapped**
******************************************************

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_StillingerWeberPairTermCapped, \
                      interaction_VerletListStillingerWeberPairTermCapped, \
                      interaction_VerletListAdressStillingerWeberPairTermCapped, \
                      interaction_VerletListHadressStillingerWeberPairTermCapped, \
                      interaction_CellListStillingerWeberPairTermCapped, \
                      interaction_FixedPairListStillingerWeberPairTermCapped

class StillingerWeberPairTermCappedLocal(PotentialLocal, interaction_StillingerWeberPairTermCapped):
  'The (local) Lennard-Jones potential.'
  def __init__(self, A, B, p, q, epsilon=1.0, sigma=1.0, cutoff=infinity, caprad = 0.0):
    """Initialize the local Lennard Jones object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_StillingerWeberPairTermCapped, A, B, p, q, epsilon, sigma, cutoff, caprad)

class VerletListStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_VerletListStillingerWeberPairTermCapped):
  'The (local) Lennard Jones interaction using Verlet lists.'
  def __init__(self, vl):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListStillingerWeberPairTermCapped, vl)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

  def getPotential(self, type1, type2):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getPotential(self, type1, type2)

  def getVerletListLocal(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getVerletList(self)
    
  def getCaprad(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getCaprad(self)

class VerletListAdressStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_VerletListAdressStillingerWeberPairTermCapped):
  'The (local) Lennard Jones interaction using Verlet lists.'
  def __init__(self, vl, fixedtupleList):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListAdressStillingerWeberPairTermCapped, vl, fixedtupleList)

  def setPotentialAT(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialAT(self, type1, type2, potential)

  def setPotentialCG(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialCG(self, type1, type2, potential)
      
class VerletListHadressStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_VerletListHadressStillingerWeberPairTermCapped):
  'The (local) Lennard Jones interaction using Verlet lists.'
  def __init__(self, vl, fixedtupleList):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListHadressStillingerWeberPairTermCapped, vl, fixedtupleList)

  def setPotentialAT(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialAT(self, type1, type2, potential)

  def setPotentialCG(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialCG(self, type1, type2, potential)

class CellListStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_CellListStillingerWeberPairTermCapped):
  'The (local) Lennard Jones interaction using cell lists.'
  def __init__(self, stor):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_CellListStillingerWeberPairTermCapped, stor)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListStillingerWeberPairTermCappedLocal(InteractionLocal, interaction_FixedPairListStillingerWeberPairTermCapped):
  'The (local) Lennard-Jones interaction using FixedPair lists.'
  def __init__(self, system, vl, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedPairListStillingerWeberPairTermCapped, system, vl, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class StillingerWeberPairTermCapped(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
          cls = 'espressopp.interaction.StillingerWeberPairTermCappedLocal',
          pmiproperty = ['A', 'B', 'p', 'q', 'epsilon', 'sigma', 'caprad'],
          pmiinvoke = ['getCaprad']
        )

    class VerletListStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotential', 'getPotential', 'getVerletList']
        )

    class VerletListAdressStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListAdressStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )
        
    class VerletListHadressStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListHadressStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )

    class CellListStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.CellListStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotential']
        )
        
    class FixedPairListStillingerWeberPairTermCapped(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.FixedPairListStillingerWeberPairTermCappedLocal',
          pmicall = ['setPotential']
        )
