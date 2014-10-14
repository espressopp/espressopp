#  Copyright (C) 2014
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


"""
*********************************************
**espresso.interaction.DihedralHarmonicNCos**
*********************************************

The dihedral harmonic potential

.. math::

   U(\phi_{ijkl}) = K\cdot[1+cos(N\cdot\phi_{ijkl} - \phi_0)]

where the `K` is a constant, the angles should be provided in radians.
The `N` is a multiplicity.

Reference: http://www.uark.edu/ua/fengwang/DLPOLY2/node49.html
"""


# pylint: disable=W0401, W0614, W0212
from espresso import pmi
from espresso.esutil import *

from espresso.interaction.DihedralPotential import *
from espresso.interaction.Interaction import *
# pylint: disable=F0401
from _espresso import interaction_DihedralHarmonicNCos
from _espresso import interaction_FixedQuadrupleListDihedralHarmonicNCos


class DihedralHarmonicNCosLocal(DihedralPotentialLocal, interaction_DihedralHarmonicNCos):
  'The (local) DihedralHarmoniNCos potential.'
  def __init__(self, K=0.0, phi0=0.0, multiplicity=1):
    """Initialize the local DihedralHarmonicNCos object."""
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


if pmi.isController:
  class DihedralHarmonicNCos(DihedralPotential):
    'The DihedralHarmonicNCos potential.'
    pmiproxydefs = dict(
      cls='espresso.interaction.DihedralHarmonicNCosLocal',
      pmiproperty=['K', 'phi', 'multiplicity']
    )

  class FixedQuadrupleListDihedralHarmonicNCos(Interaction):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls='espresso.interaction.FixedQuadrupleListDihedralHarmonicNCosLocal',
      pmicall=['setPotential', 'getFixedQuadrupleList']
    )
