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
Calculates the Angular Harmonic potential as:

.. math::

    U = K (\\theta - \\theta_0)^2,

where angle :math:`\\theta` is the planar angle formed by three binded particles
(triplet or triple). The usual coefficient of :math:`1/2` is included in :math:`K`.

This potential is employed by:

.. py:class:: espressopp.interaction.AngularHarmonic (K = 1.0, theta0 = 0.0)

    :param real K: energy amplitude
    :param real theta0: angle in radians
    :rtype: triple potential

A triple potential applied to every triple in the system creates an *interaction*.
This is done via:

.. py:class:: espressopp.interaction.FixedTripleListAngularHarmonic (system, fixed_triple_list, potential)

    :param shared_ptr system: system object
    :param list fixed_triple_list: a fixed list of all triples in the system
    :param potential: triple potential (in this case, :py:class:`espressopp.interaction.AngularHarmonic`).
    :rtype: interaction

    **Methods**

    .. py:method:: getFixedTripleList()

        :rtype: A Python list of fixed triples (e.g., in the chains).

    .. py:method:: setPotential(type1, type2, potential)

        :param type1:
        :param type2:
        :param potential:
        :type type1:
        :type type2:
        :type potential:

**Example 1.** Creating a fixed triple list by :py:class:`espressopp.FixedTripleList`.

    >>> # we assume a polymer solution of n_chains of the length chain_len each.
    >>> # At first, create a list_of_triples for the system:
    >>> N = n_chains * chain_len            # number of particles in the system
    >>> list_of_tripples = []               # empty list of triples
    >>> for n in range (n_chains):          # loop over chains
    >>>     for m in range (chain_len):     # loop over chain beads
    >>>         pid = n * chain_len + m
    >>>         if (pid > 1) and (pid < N - 1):
    >>>             list_of_tripples.append( (pid-1, pid, pid+1) )
    >>>
    >>> # create fixed triple list
    >>> fixed_triple_list = espressopp.FixedTripleList(system.storage)
    >>> fixed_triple_list.addTriples(list_of_triples)

**Example 2.** Employing an Angular Harmonic potential.

    >>> # Note, the fixed_triple_list has to be generated in advance! (see Example 1)
    >>>
    >>> # set up the potential
    >>> potAngHarm = espressopp.interaction.AngularHarmonic(K=0.5, theta0=0.0)
    >>>
    >>> # set up the interaction
    >>> interAngHarm = espressopp.interaction.FixedTripleListAngularHarmonic(system, fixed_triple_list, potAngHarm)
    >>>
    >>> # finally, add the interaction to the system
    >>> system.addInteraction(interAngHarm)

"""

from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.AngularPotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_AngularHarmonic, interaction_FixedTripleListAngularHarmonic

class AngularHarmonicLocal(AngularPotentialLocal, interaction_AngularHarmonic):

    def __init__(self, K=1.0, theta0=0.0):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_AngularHarmonic, K, theta0)

class FixedTripleListAngularHarmonicLocal(InteractionLocal, interaction_FixedTripleListAngularHarmonic):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedTripleListAngularHarmonic, system, vl, potential)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

if pmi.isController:
    class AngularHarmonic(AngularPotential):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.AngularHarmonicLocal',
            pmiproperty = ['K', 'theta0']
            )

    class FixedTripleListAngularHarmonic(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedTripleListAngularHarmonicLocal',
            pmicall = ['setPotential', 'getFixedTripleList']
            )
