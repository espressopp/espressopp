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
************************************************
**espresso.interaction.StillingerWeberPairTerm**
************************************************

"""
from espresso import pmi, infinity
from espresso.esutil import *

from espresso.interaction.Potential import *
from espresso.interaction.Interaction import *
from _espresso import interaction_StillingerWeberPairTerm, \
                      interaction_VerletListStillingerWeberPairTerm, \
                      interaction_VerletListAdressStillingerWeberPairTerm, \
                      interaction_VerletListHadressStillingerWeberPairTerm, \
                      interaction_CellListStillingerWeberPairTerm, \
                      interaction_FixedPairListStillingerWeberPairTerm

class StillingerWeberPairTermLocal(PotentialLocal, interaction_StillingerWeberPairTerm):
  'The (local) Lennard-Jones potential.'
  def __init__(self, A, B, p, q, epsilon=1.0, sigma=1.0, cutoff=infinity):
    """Initialize the local Lennard Jones object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_StillingerWeberPairTerm, A, B, p, q, epsilon, sigma, cutoff)

class VerletListStillingerWeberPairTermLocal(InteractionLocal, interaction_VerletListStillingerWeberPairTerm):
  'The (local) Lennard Jones interaction using Verlet lists.'
  def __init__(self, vl):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListStillingerWeberPairTerm, vl)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

  def getPotential(self, type1, type2):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getPotential(self, type1, type2)

  def getVerletListLocal(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getVerletList(self)

class VerletListAdressStillingerWeberPairTermLocal(InteractionLocal, interaction_VerletListAdressStillingerWeberPairTerm):
  'The (local) Lennard Jones interaction using Verlet lists.'
  def __init__(self, vl, fixedtupleList):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListAdressStillingerWeberPairTerm, vl, fixedtupleList)

  def setPotentialAT(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialAT(self, type1, type2, potential)

  def setPotentialCG(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialCG(self, type1, type2, potential)
      
class VerletListHadressStillingerWeberPairTermLocal(InteractionLocal, interaction_VerletListHadressStillingerWeberPairTerm):
  'The (local) Lennard Jones interaction using Verlet lists.'
  def __init__(self, vl, fixedtupleList, KTI = False):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListHadressStillingerWeberPairTerm, vl, fixedtupleList, KTI)

  def setPotentialAT(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialAT(self, type1, type2, potential)

  def setPotentialCG(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotentialCG(self, type1, type2, potential)
      
class CellListStillingerWeberPairTermLocal(InteractionLocal, interaction_CellListStillingerWeberPairTerm):
  'The (local) Lennard Jones interaction using cell lists.'
  def __init__(self, stor):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_CellListStillingerWeberPairTerm, stor)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListStillingerWeberPairTermLocal(InteractionLocal, interaction_FixedPairListStillingerWeberPairTerm):
  'The (local) Lennard-Jones interaction using FixedPair lists.'
  def __init__(self, system, vl, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedPairListStillingerWeberPairTerm, system, vl, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class StillingerWeberPairTerm(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
          cls = 'espresso.interaction.StillingerWeberPairTermLocal',
          pmiproperty = ['A', 'B', 'p', 'q', 'epsilon', 'sigma']
        )

    class VerletListStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.VerletListStillingerWeberPairTermLocal',
          pmicall = ['setPotential', 'getPotential', 'getVerletList']
        )

    class VerletListAdressStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.VerletListAdressStillingerWeberPairTermLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )
        
    class VerletListHadressStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.VerletListHadressStillingerWeberPairTermLocal',
          pmicall = ['setPotentialAT', 'setPotentialCG']
        )

    class CellListStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.CellListStillingerWeberPairTermLocal',
          pmicall = ['setPotential']
        )
        
    class FixedPairListStillingerWeberPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espresso.interaction.FixedPairListStillingerWeberPairTermLocal',
          pmicall = ['setPotential']
        )
