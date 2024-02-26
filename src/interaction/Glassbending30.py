#  Copyright (C) 2012,2013,2016,2018
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
Calculates the Glassbending30 potential as:

Reference: Hsiao-Ping Hsu and Kurt Kremer
           "A coarse-grained polymer model for studying the glass transition",
           J. Chem. Phys. 150, 091101 (2019)
.. math::

   U = -Ba \sin^2(Bb*(pi-\theta) = -Ba/2.*(1-\cos (2*Bb*(pi-\theta)))

.. function:: espressopp.interaction.Glassbending30(K, theta0)

      :param Ba: (default: 0.0)
      :param Bb: (default: 0.0)
      :param Bbthreshold: (default: 0.0)
      :type Ba: real
      :type Bb: real
      :type Bbthreshold: real


.. function:: espressopp.interaction.FixedTripleListGlassbending30(system, vl, potential)

                :param system:
                :param vl:
                :param potential:
                :type system:
                :type vl:
                :type potential:

.. function:: espressopp.interaction.FixedTripleListGlassbending30.getFixedTripleList()

                :rtype: A Python list of lists.

.. function:: espressopp.interaction.FixedTripleListGlassbending30.setPotential(potential)

                :param potential:
                :type potential:
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.AngularPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_Glassbending30, interaction_FixedTripleListGlassbending30

class Glassbending30Local(AngularPotentialLocal, interaction_Glassbending30):

    def __init__(self, Ba=1.0, Bb=0.0, Bbthreshold=0.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_Glassbending30, Ba, Bb, Bbthreshold)

class FixedTripleListGlassbending30Local(InteractionLocal, interaction_FixedTripleListGlassbending30):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListGlassbending30, system, vl, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

    def getFixedTripleList(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getFixedTripleList(self)

if pmi.isController:
    class Glassbending30(AngularPotential):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.Glassbending30Local',
            pmiproperty = ['Ba', 'Bb', 'Bbthreshold']
            )

    class FixedTripleListGlassbending30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListGlassbending30Local',
            pmicall = ['setPotential','getFixedTripleList']
            )
