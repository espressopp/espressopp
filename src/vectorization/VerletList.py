#  Copyright (C) 2019
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
*********************
espressopp.vectorization.VerletList
*********************

LIMITATION: This class does not support exclusion lists as they are not supported in the
vectorized version of VerletList.

.. function:: espressopp.vectorization.VerletList(system, vec, cutoff, exclusionlist, build_order)

		:param system:
		:param vec: Vectorization object
		:param cutoff:
		:param exclusionlist: (default: [])
		:param build_order: (default: 0) Whether to build the verlet list in the p->nc or nc->p order
		:type system:
		:type vec:
		:type cutoff:
		:type exclusionlist:
		:type build_order:

.. function:: espressopp.vectorization.VerletList.exclude(exclusionlist)

		:param exclusionlist:
		:type exclusionlist:
		:rtype:

		Throws an error since exclusion lists are not supported.

.. function:: espressopp.vectorization.VerletList.getAllPairs()

		:rtype:

.. function:: espressopp.vectorization.VerletList.localSize()

		:rtype:

.. function:: espressopp.vectorization.VerletList.totalSize()

		:rtype:
"""
from espressopp import pmi
import _espressopp
import espressopp
from espressopp.esutil import cxxinit

class VerletListLocal(_espressopp.vectorization_VerletList):


    def __init__(self, system, vec, cutoff, exclusionlist=[], build_order=0):

        if pmi.workerIsActive():
            if (exclusionlist == []):
                # rebuild list in constructor
                cxxinit(self, _espressopp.vectorization_VerletList, system, vec, cutoff, True, build_order)
            else:
                # do not rebuild list in constructor
                cxxinit(self, _espressopp.vectorization_VerletList, system, vec, cutoff, False, build_order)
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

    def exclude(self, exclusionlist):
        """
        Each processor takes the broadcasted exclusion list
        and adds it to its list.
        """
        if pmi.workerIsActive():
            for pair in exclusionlist:
                pid1, pid2 = pair
                self.cxxclass.exclude(self, pid1, pid2)
            # rebuild list with exclusions
            self.cxxclass.rebuild(self)

    def getAllPairs(self):

        if pmi.workerIsActive():
            pairs=[]
            npairs=self.localSize()
            for i in xrange(npairs):
              pair=self.cxxclass.getPair(self, i+1)
              pairs.append(pair)
            return pairs


if pmi.isController:
  class VerletList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls = 'espressopp.vectorization.VerletListLocal',
      pmiproperty = [ 'builds' ],
      pmicall = [ 'totalSize', 'exclude', 'connect', 'disconnect', 'getVerletCutoff', 'resetTimers','rebuildPairs'],
      pmiinvoke = [ 'getAllPairs','getTimers' ]
    )
