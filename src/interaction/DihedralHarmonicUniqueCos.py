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
**************************************************
**espressopp.interaction.DihedralHarmonicUniqueCos**
**************************************************

"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.DihedralUniquePotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_DihedralHarmonicUniqueCos, \
                      interaction_FixedQuadrupleAngleListDihedralHarmonicUniqueCos

class DihedralHarmonicUniqueCosLocal(DihedralUniquePotentialLocal, interaction_DihedralHarmonicUniqueCos):
  'The (local) DihedralHarmonicUniqueCos potential.'
  def __init__(self, K=0.0):
    """Initialize the local DihedralHarmonicUniqueCos object."""
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_DihedralHarmonicUniqueCos, K)

class FixedQuadrupleAngleListDihedralHarmonicUniqueCosLocal(InteractionLocal, interaction_FixedQuadrupleAngleListDihedralHarmonicUniqueCos):
  'The (local) DihedralHarmonicUniqueCos interaction using FixedQuadruple lists.'
  def __init__(self, system, fqal, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedQuadrupleAngleListDihedralHarmonicUniqueCos, system, fqal, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)
      
  def getFixedQuadrupleList(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getFixedQuadrupleAngleList(self)


if pmi.isController:
  class DihedralHarmonicUniqueCos(DihedralUniquePotential):
    'The DihedralHarmonicUniqueCos potential.'
    pmiproxydefs = dict(
      cls = 'espressopp.interaction.DihedralHarmonicUniqueCosLocal',
      pmiproperty = ['K']
    )

  class FixedQuadrupleAngleListDihedralHarmonicUniqueCos(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.interaction.FixedQuadrupleAngleListDihedralHarmonicUniqueCosLocal',
      pmicall = ['setPotential', 'getFixedQuadrupleList']
    )
