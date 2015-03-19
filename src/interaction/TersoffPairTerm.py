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
****************************************
**espressopp.interaction.TersoffPairTerm**
****************************************

"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_TersoffPairTerm, \
                      interaction_VerletListTersoffPairTerm, \
                      interaction_CellListTersoffPairTerm, \
                      interaction_FixedPairListTersoffPairTerm

class TersoffPairTermLocal(PotentialLocal, interaction_TersoffPairTerm):
  'The (local) Lennard-Jones potential.'
  def __init__(self, A, lambda1, R, D, cutoff=infinity):
    """Initialize the local Lennard Jones object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_TersoffPairTerm, A, lambda1, R, D, cutoff)

class VerletListTersoffPairTermLocal(InteractionLocal, interaction_VerletListTersoffPairTerm):
  'The (local) Lennard Jones interaction using Verlet lists.'
  def __init__(self, vl):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_VerletListTersoffPairTerm, vl)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

  def getPotential(self, type1, type2):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getPotential(self, type1, type2)

  def getVerletListLocal(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getVerletList(self)

class CellListTersoffPairTermLocal(InteractionLocal, interaction_CellListTersoffPairTerm):
  'The (local) Lennard Jones interaction using cell lists.'
  def __init__(self, stor):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_CellListTersoffPairTerm, stor)

  def setPotential(self, type1, type2, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListTersoffPairTermLocal(InteractionLocal, interaction_FixedPairListTersoffPairTerm):
  'The (local) Lennard-Jones interaction using FixedPair lists.'
  def __init__(self, system, vl, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedPairListTersoffPairTerm, system, vl, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class TersoffPairTerm(Potential):
        'The Lennard-Jones potential.'
        pmiproxydefs = dict(
          cls = 'espressopp.interaction.TersoffPairTermLocal',
          pmiproperty = ['A', 'lambda1', 'R', 'D']
        )

    class VerletListTersoffPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.VerletListTersoffPairTermLocal',
          pmicall = ['setPotential', 'getPotential', 'getVerletList']
        )

    class CellListTersoffPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.CellListTersoffPairTermLocal',
          pmicall = ['setPotential']
        )
        
    class FixedPairListTersoffPairTerm(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.interaction.FixedPairListTersoffPairTermLocal',
          pmicall = ['setPotential']
        )
