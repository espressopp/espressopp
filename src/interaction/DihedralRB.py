#  Copyright (C) 2015
#      Jakub Krajniak (jkrajniak at gmail.com)
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
*****************************
**espressopp.interaction.DihedralRB**
*****************************

The proper dihedral with Ryckaert-Bellemans form.

.. math::

    U_{rb}(\phi_{ijkl}) = \sum_{n=0}^{5} K_n (cos(\\theta))^n

where the :math:`\\theta = \phi - 180^\circ` and :math:`K_{0\dotso5}` are the coefficients.

By default the IUPAC confention is used, where :math:`\phi` is the angle between planes
:math:`ijk` and :math:`jkl`. The :math:`0^\circ` corresponds to the *cis* configuration.

Reference: http://www.gromacs.org/Documentation/Manual
"""

from espressopp import pmi
from espressopp.esutil import *  # NOQA

from espressopp.interaction.DihedralPotential import *  # NOQA
from espressopp.interaction.Interaction import *  # NOQA
from _espressopp import interaction_DihedralRB
from _espressopp import interaction_FixedQuadrupleListDihedralRB


class DihedralRBLocal(DihedralPotentialLocal, interaction_DihedralRB):
    """The (local) DihedralRB potential."""
    def __init__(self, K0=0.0, K1=0.0, K2=0.0, K3=0.0, K4=0.0, K5=0.0, iupac=True):
        """Initialize the local DihedralRB object."""
        if (not (pmi._PMIComm and pmi._PMIComm.isActive())
                or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            if iupac:
                cxxinit(self, interaction_DihedralRB, K0, -1*K1, K2, -1*K3, K4, -1*K5)
            else:
                cxxinit(self, interaction_DihedralRB, K0, K1, K2, K3, K4, K5)


class FixedQuadrupleListDihedralRBLocal(InteractionLocal, interaction_FixedQuadrupleListDihedralRB):
    """The (local) DihedralRB interaction using FixedQuadruple lists."""
    def __init__(self, system, vl, potential):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive())
                or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            cxxinit(self, interaction_FixedQuadrupleListDihedralRB, system, vl, potential)

    def setPotential(self, type1, type2, potential):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive())
                or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getFixedQuadrupleList(self):
        if (not (pmi._PMIComm and pmi._PMIComm.isActive())
                or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup()):
            return self.cxxclass.getFixedQuadrupleList(self)

if pmi.isController:
    class DihedralRB(DihedralPotential):
        """The  Ryckaert-Bellemans function for proper dihedrals.

        Args:
            K0, K1, K2, K3, K4, K5: The parameters for potential.
            iupac: If set to true then IUPAC convention for dihedrals is used (by default).
        """
        pmiproxydefs = dict(
            cls='espressopp.interaction.DihedralRBLocal'
            )

    class FixedQuadrupleListDihedralRB(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.interaction.FixedQuadrupleListDihedralRBLocal',
            pmicall=['setPotential', 'getFixedQuadrupleList']
            )
