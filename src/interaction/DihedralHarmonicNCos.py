#  Copyright (C) 2014,2016
#      Jakub Krajniak
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
*******************************************
espressopp.interaction.DihedralHarmonicNCos
*******************************************

The dihedral harmonic potential

.. math::

   U(\phi_{ijkl}) = K[1+cos(N\cdot\phi_{ijkl} - \phi_0)]

where the `K` is a constant, the angles should be provided in radians.
The `N` is a multiplicity.

Reference: http://www.uark.edu/ua/fengwang/DLPOLY2/node49.html






.. function:: espressopp.interaction.DihedralHarmonicNCos(K, phi0, multiplicity)

		:param K: (default: 0.0)
		:param phi0: (default: 0.0)
		:param multiplicity: (default: 1)
		:type K: real
		:type phi0: real
		:type multiplicity: int

.. function:: espressopp.interaction.FixedQuadrupleListDihedralHarmonicNCos(system, fql, potential)

		:param system: 
		:param fql: 
		:param potential: 
		:type system: 
		:type fql: 
		:type potential: 

.. function:: espressopp.interaction.FixedQuadrupleListDihedralHarmonicNCos.getFixedQuadrupleList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedQuadrupleListDihedralHarmonicNCos.setPotential(potential)

		:param potential: 
		:type potential: 
"""


# pylint: disable=W0401, W0614, W0212
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.DihedralPotential import *
from espressopp.interaction.Interaction import *
# pylint: disable=F0401
from _espressopp import interaction_DihedralHarmonicNCos
from _espressopp import interaction_FixedQuadrupleListDihedralHarmonicNCos
from _espressopp import interaction_FixedQuadrupleListTypesDihedralHarmonicNCos


class DihedralHarmonicNCosLocal(DihedralPotentialLocal, interaction_DihedralHarmonicNCos):

  def __init__(self, K=0.0, phi0=0.0, multiplicity=1):

    # pylint: disable=W0212
    if (not (pmi._PMIComm and pmi._PMIComm.isActive())
        or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
      cxxinit(self, interaction_DihedralHarmonicNCos, K, phi0, multiplicity)


class FixedQuadrupleListDihedralHarmonicNCosLocal(
    InteractionLocal,
    interaction_FixedQuadrupleListDihedralHarmonicNCos):
  'The (local) DihedralHarmonicNCos interaction using FixedQuadruple lists.'
  def __init__(self, system, fql, potential):
    if (not (pmi._PMIComm and pmi._PMIComm.isActive())
        or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
      cxxinit(self, interaction_FixedQuadrupleListDihedralHarmonicNCos, system, fql, potential)

  def setPotential(self, potential):
    if (not (pmi._PMIComm and pmi._PMIComm.isActive())
        or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
      self.cxxclass.setPotential(self, potential)

  def getFixedQuadrupleList(self):
    if (not (pmi._PMIComm and pmi._PMIComm.isActive())
        or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
      return self.cxxclass.getFixedQuadrupleList(self)


class FixedQuadrupleListTypesDihedralHarmonicNCosLocal(InteractionLocal,
                                                   interaction_FixedQuadrupleListTypesDihedralHarmonicNCos):
  def __init__(self, system, fql):
    if pmi.workerIsActive():
      cxxinit(self, interaction_FixedQuadrupleListTypesDihedralHarmonicNCos, system, fql)

  def setPotential(self, type1, type2, type3, type4, potential):
    if pmi.workerIsActive():
      self.cxxclass.setPotential(self, type1, type2, type3, type4, potential)

  def getPotential(self, type1, type2, type3, type4):
    if pmi.workerIsActive():
      return self.cxxclass.getPotential(self, type1, type2, type3, type4)


if pmi.isController:
  class DihedralHarmonicNCos(DihedralPotential):
    'The DihedralHarmonicNCos potential.'
    pmiproxydefs = dict(
      cls='espressopp.interaction.DihedralHarmonicNCosLocal',
      pmiproperty=['K', 'phi0', 'multiplicity']
    )

  class FixedQuadrupleListDihedralHarmonicNCos(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls='espressopp.interaction.FixedQuadrupleListDihedralHarmonicNCosLocal',
      pmicall=['setPotential', 'getFixedQuadrupleList']
    )

  class FixedQuadrupleListTypesDihedralHarmonicNCos(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.interaction.FixedQuadrupleListTypesDihedralHarmonicNCosLocal',
      pmicall = ['setPotential','getPotential','setFixedQuadrupleList','getFixedQuadrupleList']
    )
