#  Copyright (C) 2014
#      Pierre de Buyl
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
*******************************************
**espressopp.interaction.LennardJones93Wall**
*******************************************

This class defines a Lennard-Jones 9-3 SingleParticlePotential in the direction x.

.. math:: V(r) = \\epsilon \\left( \\left(\\frac{\sigma}{r}\\right)^9 - \\left(\\frac{\sigma}{r}\\right)^3 \\right)

where :math:`r` is the distance from the lower or upper wall in the x
direction. :math:`V(r)=0` after a distance `sigmaCutoff`.

The parameters have to be defined for every species present in the system with
`setParams` and can be retrieved with `getParams`.

Example:

    >>> LJ93 = espressopp.interaction.LennardJones93Wall()
    >>> LJ93.setParams(0, 6., 1., wall_cutoff)
    >>> SPLJ93 = espressopp.interaction.SingleParticleLennardJones93Wall(system, LJ93)
    >>> system.addInteraction(SPLJ93)

"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.SingleParticlePotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJones93Wall, interaction_SingleParticleLennardJones93Wall


class LennardJones93WallLocal(SingleParticlePotentialLocal, interaction_LennardJones93Wall):
    'The (local) LennardJones93Wall potential.'
    def __init__(self):
        """Initialize the local LennardJones93Wall object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_LennardJones93Wall)
    def getParams(self, type_var):
        """Return the epsilon, sigma, sigmaCutoff and r0 parameters for particles of type `type_var`.
        """
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getParams(self, type_var)

    def setParams(self, type_var, epsilon, sigma, sigmaCutoff, r0):
        """Set the epsilon, sigma, sigmaCutoff and r0 parameters for particles of type `type_var`.
        """
        self.cxxclass.setParams(self, type_var, epsilon, sigma, sigmaCutoff, r0)

class SingleParticleLennardJones93WallLocal(InteractionLocal, interaction_SingleParticleLennardJones93Wall):
    'The (local) LennardJones93Wall interaction.'
    def __init__(self, system, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_SingleParticleLennardJones93Wall, system, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class LennardJones93Wall(SingleParticlePotential):
        'The LennardJones93Wall potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJones93WallLocal',
            pmicall = ['setParams', 'getParams']
            )

    class SingleParticleLennardJones93Wall(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.SingleParticleLennardJones93WallLocal',
            pmicall = ['setPotential']
            )
