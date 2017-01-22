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
espressopp.interaction.DihedralHarmonicCos
******************************************

.. math::

	U = K (cos(\phi) - cos(\phi_0))^2

.. function:: espressopp.interaction.DihedralHarmonicCos(K, phi0)

		:param K: (default: 0.0)
		:param phi0: (default: 0.0)
		:type K: real
		:type phi0: real

.. function:: espressopp.interaction.FixedQuadrupleListDihedralHarmonicCos(system, fql, potential)

		:param system: 
		:param fql: 
		:param potential: 
		:type system: 
		:type fql: 
		:type potential: 

.. function:: espressopp.interaction.FixedQuadrupleListDihedralHarmonicCos.getFixedQuadrupleList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedQuadrupleListDihedralHarmonicCos.setPotential(potential)

		:param potential: 
		:type potential: 
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.DihedralPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_DihedralHarmonicCos, interaction_FixedQuadrupleListDihedralHarmonicCos

class DihedralHarmonicCosLocal(DihedralPotentialLocal, interaction_DihedralHarmonicCos):

  def __init__(self, K=0.0, phi0=0.0):

    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_DihedralHarmonicCos, K, phi0)

class FixedQuadrupleListDihedralHarmonicCosLocal(InteractionLocal, interaction_FixedQuadrupleListDihedralHarmonicCos):

  def __init__(self, system, fql, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, interaction_FixedQuadrupleListDihedralHarmonicCos, system, fql, potential)

  def setPotential(self, potential):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.setPotential(self, potential)
      
  def getFixedQuadrupleList(self):
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      return self.cxxclass.getFixedQuadrupleList(self)


if pmi.isController:
  class DihedralHarmonicCos(DihedralPotential):
    'The DihedralHarmonicCos potential.'
    pmiproxydefs = dict(
      cls = 'espressopp.interaction.DihedralHarmonicCosLocal',
      pmiproperty = ['K', 'phi']
    )

  class FixedQuadrupleListDihedralHarmonicCos(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.interaction.FixedQuadrupleListDihedralHarmonicCosLocal',
      pmicall = ['setPotential', 'getFixedQuadrupleList']
    )
