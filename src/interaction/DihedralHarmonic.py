#  Copyright (C) 2012,2013,2015(H),2016
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
***************************************
espressopp.interaction.DihedralHarmonic
***************************************

The dihedral harmonic potential

.. math::

   U(\phi_{ijkl}) = 0.5K[\phi_{ijkl} - \phi_0)]^2

where the `K` is a constant, the angles should be provided in radians.

Reference: Gromacs Manual 4.6.1, section 4.2.11 (page 79-80), equation 4.60



.. function:: espressopp.interaction.DihedralHarmonic(K, phi0)

		:param K: (default: 0.0)
		:param phi0: (default: 0.0)
		:type K: real
		:type phi0: real

.. function:: espressopp.interaction.FixedQuadrupleListDihedralHarmonic(system, fql, potential)

		:param system: 
		:param fql: 
		:param potential: 
		:type system: 
		:type fql: 
		:type potential: 

.. function:: espressopp.interaction.FixedQuadrupleListDihedralHarmonic.getFixedQuadrupleList()

		:rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedQuadrupleListDihedralHarmonic.setPotential(potential)

		:param potential: 
		:type potential: 

**Example of usage**

>>> # The following example shows how to add a torsional potential to particles 1,2,3,4
>>> fql = espressopp.FixedQuadrupleList(system.storage)
>>> fql.addQuadruples([(1,2,3,4)])
>>> #phi0 is in radians, IUPAC convention definition
>>> interaction = espressopp.interaction.FixedQuadrupleListDihedralHarmonic(system,fql,potential=espressopp.interaction.DihedralHarmonic(K=1.0,phi0=0.0))
>>> system.addInteraction(interaction)

"""


# pylint: disable=W0401, W0614, W0212
from espressopp.esutil import *

from espressopp.interaction.DihedralPotential import *
from espressopp.interaction.Interaction import *
# pylint: disable=F0401
from _espressopp import interaction_DihedralHarmonic
from _espressopp import interaction_FixedQuadrupleListDihedralHarmonic
from _espressopp import interaction_FixedQuadrupleListTypesDihedralHarmonic


class DihedralHarmonicLocal(DihedralPotentialLocal, interaction_DihedralHarmonic):

  def __init__(self, K=0.0, phi0=0.0):

    # pylint: disable=W0212
    if (not (pmi._PMIComm and pmi._PMIComm.isActive())
        or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
      cxxinit(self, interaction_DihedralHarmonic, K, phi0)


class FixedQuadrupleListDihedralHarmonicLocal(
    InteractionLocal,
    interaction_FixedQuadrupleListDihedralHarmonic):
  'The (local) DihedralHarmonic interaction using FixedQuadruple lists.'
  def __init__(self, system, fql, potential):
    if (not (pmi._PMIComm and pmi._PMIComm.isActive())
        or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
      cxxinit(self, interaction_FixedQuadrupleListDihedralHarmonic, system, fql, potential)

  def setPotential(self, potential):
    if (not (pmi._PMIComm and pmi._PMIComm.isActive())
        or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
      self.cxxclass.setPotential(self, potential)

  def getFixedQuadrupleList(self):
    if (not (pmi._PMIComm and pmi._PMIComm.isActive())
        or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
      return self.cxxclass.getFixedQuadrupleList(self)


class FixedQuadrupleListTypesDihedralHarmonicLocal(InteractionLocal,
                                                    interaction_FixedQuadrupleListTypesDihedralHarmonic):
  def __init__(self, system, fql):
    if pmi.workerIsActive():
      cxxinit(self, interaction_FixedQuadrupleListTypesDihedralHarmonic, system, fql)

  def setPotential(self, type1, type2, type3, type4, potential):
    if pmi.workerIsActive():
      self.cxxclass.setPotential(self, type1, type2, type3, type4, potential)

  def getPotential(self, type1, type2, type3, type4):
    if pmi.workerIsActive():
      return self.cxxclass.getPotential(self, type1, type2, type3, type4)


if pmi.isController:
  class DihedralHarmonic(DihedralPotential):
    'The DihedralHarmonic potential.'
    pmiproxydefs = dict(
      cls='espressopp.interaction.DihedralHarmonicLocal',
      pmiproperty=['K', 'phi0']
    )

  class FixedQuadrupleListDihedralHarmonic(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls='espressopp.interaction.FixedQuadrupleListDihedralHarmonicLocal',
      pmicall=['setPotential', 'getFixedQuadrupleList']
    )

  class FixedQuadrupleListTypesDihedralHarmonic(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.interaction.FixedQuadrupleListTypesDihedralHarmonicLocal',
      pmicall = ['setPotential','getPotential','setFixedQuadrupleList','getFixedQuadrupleList']
    )
