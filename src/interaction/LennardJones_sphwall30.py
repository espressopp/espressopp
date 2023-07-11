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

r"""
*****************************************
espressopp.interaction.LennardJones_sphwall30
*****************************************

This class defines a Lennard-Jones SingleParticlePotential in spherical confinement.

.. math:: V(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{R-r_s} \right)^{12} -
        \left( \frac{\sigma}{R-r_s} \right)^{6}

where :math:`R` is the radius of sphere, r_s is the distance between the monomer and center of sphere
direction. :math:`V(r)=0` after a distance `r_{cut}`.

The parameters have to be defined for every species present in the system with
`setParams` and can be retrieved with `getParams`.

    >>> LJ_sphw = espressopp.interaction.LennardJones_sphwall30()
    >>> LJ_sphw.setParams(0, eps_wall, sigma, radius, lcutoff)   
    >>> SPLJ_sphw = espressopp.interaction.SingleParticleLennardJones_sphwall30(system, LJ_sphw)
    >>> system.addInteraction(SPLJ_sphw)


.. function:: espressopp.interaction.LennardJones_sphwall30()


.. function:: espressopp.interaction.LennardJones_sphwall30.getParams(type_var)

                :param type_var:
                :type type_var:
                :rtype:

.. function:: espressopp.interaction.nLennardJones_sphwall30.setParams(type_var, epsilon, sigma, radius, lcutoff)

                :param type_var:
                :param epsilon:
                :param sigma:
                :param radius:
                :param lcutoff:
                :type type_var:
                :type epsilon:
                :type sigma:
                :type radius:
                :type lcutoff:


.. function:: espressopp.interaction.SingleParticleLennardJones_sphwall30(system, potential)

                :param system:
                :param potential:
                :type system:
                :type potential:

.. function:: espressopp.interaction.SingleParticleLennardJones_sphwall30.setPotential(potential)

                :param potential:
                :type potential:
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.SingleParticlePotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJones_sphwall30, interaction_SingleParticleLennardJones_sphwall30


class LennardJones_sphwall30Local(SingleParticlePotentialLocal, interaction_LennardJones_sphwall30):

    def __init__(self):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_LennardJones_sphwall30)
    def getParams(self, type_var):


        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getParams(self, type_var)

    def setParams(self, type_var, epsilon, sigma, radius, lcutoff):


        self.cxxclass.setParams(self, type_var, epsilon, sigma, radius, lcutoff)

class SingleParticleLennardJones_sphwall30Local(InteractionLocal, interaction_SingleParticleLennardJones_sphwall30):

    def __init__(self, system, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_SingleParticleLennardJones_sphwall30, system, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class LennardJones_sphwall30(SingleParticlePotential):
        'The LennardJones_sphwall30 potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJones_sphwall30Local',
            pmicall = ['setParams', 'getParams']
            )

    class SingleParticleLennardJones_sphwall30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.SingleParticleLennardJones_sphwall30Local',
            pmicall = ['setPotential']
            )
