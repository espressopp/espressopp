#  Copyright (C) 2019-2022
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
espressopp.hpx4espp.VerletList
***********************************

NOTE:
    - This vectorized version of Verlet List does not support exclusion lists.
    - Only the force calculation is vectorized. Calculating the energy and virial still rely
      on the original Particle pair list so rebuildPairs() needs to be called before any analysis.

.. function:: espressopp.hpx4espp.VerletList(system, vec, cutoff, exclusionlist, build_order)

		:param system:
		:param vec: Vectorization object
		:param cutoff:
		:param exclusionlist: (default: [])
		:type system:
		:type vec:
		:type cutoff:
		:type exclusionlist:

.. function:: espressopp.hpx4espp.VerletList.exclude(exclusionlist)

		:param exclusionlist:
		:type exclusionlist:
		:rtype:

		Throws an error since exclusion lists are not supported.

.. function:: espressopp.hpx4espp.VerletList.getAllPairs()

		:rtype:

.. function:: espressopp.hpx4espp.VerletList.localSize()

		:rtype:

.. function:: espressopp.hpx4espp.VerletList.totalSize()

		:rtype:

.. function:: espressopp.hpx4espp.VerletList.rebuildPairs()

		Rebuilds the non-vectorized Particle pair list which is needed when calculating the energy
		and virialduring the analysis stages.

"""
from espressopp import pmi
from _espressopp import hpx4espp_VerletList
import espressopp
from espressopp.esutil import cxxinit

class VerletListLocal(hpx4espp_VerletList):


    def __init__(self, system, cutoff, exclusionlist=[]):

        if pmi.workerIsActive():
            resortOnLoad = not system.storage.connectedOffload()
            if (exclusionlist == []):
                # rebuild list in constructor
                cxxinit(self, hpx4espp_VerletList, system, system.storage, cutoff, True, resortOnLoad)
            else:
                # do not rebuild list in constructor
                cxxinit(self, hpx4espp_VerletList, system, system.storage, cutoff, False, resortOnLoad)
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
            cls = 'espressopp.hpx4espp.VerletListLocal',
            pmiproperty = [ 'builds' ],
            pmicall = [ 'totalSize', 'exclude', 'connect', 'disconnect', 'getVerletCutoff', 'resetTimers','rebuildPairs', 'preallocFactor'],
            pmiinvoke = [ 'getTimers', 'localSize', 'getSubSizes']
            #   pmiinvoke = [ 'getAllPairs','getTimers', 'localSize' ]
        )
