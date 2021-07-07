#  Copyright (C) 2019-2021
#      Max Planck Institute for Polymer Research & JGU Mainz
#  Copyright (C) 2012,2013,2015,2016
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
***********************************
espressopp.vec.VerletList
***********************************

NOTE:
    - This vectorized version of Verlet List does not support exclusion lists.
    - Only the force calculation is vectorized. Calculating the energy and virial still rely
      on the original Particle pair list so rebuildPairs() needs to be called before any analysis.

.. function:: espressopp.vec.VerletList(system, vec, cutoff, exclusionlist, build_order)

		:param system:
		:param vec: Vectorization object
		:param cutoff:
		:param exclusionlist: (default: [])
		:type system:
		:type vec:
		:type cutoff:
		:type exclusionlist:

.. function:: espressopp.vec.VerletList.exclude(exclusionlist)

		:param exclusionlist:
		:type exclusionlist:
		:rtype:

		Throws an error since exclusion lists are not supported.

.. function:: espressopp.vec.VerletList.getAllPairs()

		:rtype:

.. function:: espressopp.vec.VerletList.localSize()

		:rtype:

.. function:: espressopp.vec.VerletList.totalSize()

		:rtype:

.. function:: espressopp.vec.VerletList.rebuildPairs()

		Rebuilds the non-vectorized Particle pair list which is needed when calculating the energy
		and virialduring the analysis stages.

"""
from espressopp import pmi
from _espressopp import vec_VerletList
import espressopp
from espressopp.esutil import cxxinit

class VerletListLocal(vec_VerletList):

    def __init__(self, system, cutoff, exclusionlist=[]):
        if pmi.workerIsActive():

            if (exclusionlist == []):
                # rebuild list in constructor
                cxxinit(self, vec_VerletList, system, cutoff, True)
            else:
                # do not rebuild list in constructor
                cxxinit(self, vec_VerletList, system, cutoff, False)
                # add exclusions
                for pair in exclusionlist:
                    pid1, pid2 = pair
                    self.cxxclass.exclude(self, pid1, pid2)
                # now rebuild list with exclusions
                self.cxxclass.rebuild(self)


    def totalSize(self):
        if pmi.workerIsActive():
            return self.cxxclass.totalSize(self)

    def localSize(self):
        if pmi.workerIsActive():
            return self.cxxclass.localSize(self)

    # def exclude(self, exclusionlist):
    #     """
    #     Each processor takes the broadcasted exclusion list
    #     and adds it to its list.
    #     """
    #     if pmi.workerIsActive():
    #         for pair in exclusionlist:
    #             pid1, pid2 = pair
    #             self.cxxclass.exclude(self, pid1, pid2)
    #         # rebuild list with exclusions
    #         self.cxxclass.rebuild(self)

    # def getAllPairs(self):

    #     if pmi.workerIsActive():
    #         pairs=[]
    #         npairs=self.localSize()
    #         for i in xrange(npairs):
    #           pair=self.cxxclass.getPair(self, i+1)
    #           pairs.append(pair)
    #         return pairs


if pmi.isController:
    class VerletList(object, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls = 'espressopp.vec.VerletListLocal',
            pmiproperty = [ 'builds' ],
            pmicall = [ 'totalSize', 'exclude', 'connect', 'disconnect', 'getVerletCutoff', 'resetTimers','rebuildPairs', 'preallocFactor'],
            pmiinvoke = [ 'getTimers', 'localSize' ]
            # pmiinvoke = [ 'getAllPairs','getTimers', 'localSize' ]
        )
