#  Copyright (C) 2017
#      Jakub Krajniak (jkrajniak at gmail.com)
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
espressopp.VerletList
*********************


.. function:: espressopp.VerletList(system, cutoff, exclusionlist)

		:param system: E++ system object
		:param cutoff: cut-off
		:param exclusionlist: exclusion list (default: None)
		:type system: espressopp.System
		:type cutoff: double
		:type exclusionlist: None, list or espressopp.DynamicExcludeList

.. function:: espressopp.VerletList.exclude(exclusionlist)

		:param exclusionlist: The static or dynamic exclusion list.
		:type exclusionlist: If static then of type `list` if dynamic then `espressopp.DynamicExcludeList`

.. function:: espressopp.VerletList.getAllPairs()

		:rtype: returns a list of all pairs

.. function:: espressopp.VerletList.localSize()

		:rtype: returns local number of pairs

.. function:: espressopp.VerletList.totalSize()

		:rtype: returns global number of pairs


******************************
espressopp.DynamicExcludeList
******************************

This special kind of exclusion list for espressopp.VerletList creates entries to exclude
by observing fixed pair lists, fixed triple lists or fixed quadruple lists.
Whenever new tuple, triplet or quadruplet is added or removed from such fixed lists, this is
reflected in the new exclusion.

In the case of fixed pair list, new exclusion is formed from the pair of new tuple.
The newly created triplet of *(pid1, pid2, pid3)* results in new exclusion *(pid1, pid3)* and
newly created quadruplet of (pid1, pid2, pid3, pid4) results in following exclusions:
*(pid1, pid3), (pid2, pid3) and (pid1, pid4)*.

.. function:: espressopp.DynamicExcludeList(integrator, exlcusionlist=None)

   :param integrator: The integrator object.
   :type integrator: espressopp.integrator.MDIntegrator
   :param exclusionlist: The optional initial list of exclusions.
   :type exclusionlist: list

.. function:: espressopp.DynamicExcludeList.exclude(pid1, pid2)

   Adds pair of `(pid1, pid2)` to the exclude list.

   :param pid1: The id of first particle.
   :type pid1: int
   :param pid2: The id of second particle.
   :type pid2: int

.. function:: espressopp.DynamicExcludeList.unexclude(pid1, pid2)

   Removes pair of `(pid1, pid2)` from the exclude list.

   :param pid1: The id of first particle.
   :type pid1: int
   :param pid2: The id of second particle.
   :type pid2: int

.. function:: espressopp.DynamicExcludeList.observe_tuple(fpl)

   Add `espressopp.FixedPairList` to the observation list. Whenever new pair
   will be added to that object, then this pair will also be added to the excluded list.

   :param fpl: The fixed pair list object.
   :type fpl: espressopp.FixedPairList

.. function:: espressopp.DynamicExcludeList.observe_triple(ftl)

   It is similar to the `observe_tuple` command, but observes `espressopp.FixedTripleList`.
   Whenever new triplet is added `(pid1, pid2, pid3)` then the pairs `(pid1, pid2), (pid2, pid3), (pid1, pid3)`
   are added to exclude list.

   :param ftl: The fixed triple list object.
   :type ftl: espressopp.FixedTripleList

.. function:: espressopp.DynamicExcludeList.observe_quadruple(fql)

   Observe `espressopp.FixedQuadrupleList`. The pairs `(pid1, pid2), (pid1, pid3), (pid1, pid4), (pid2, pid3)`,
   `(pid2, pid4), (pid3, pid4)` are added to exclude list.

   :param fql: The fixed quadruple list object.
   :type fql: espressopp.FixedQuadrupleList

"""
from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit


class DynamicExcludeListLocal(_espressopp.DynamicExcludeList):
    def __init__(self, integrator, exclusionlist=None):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.DynamicExcludeList, integrator)
            if exclusionlist is not None:
                for pid1, pid2 in exclusionlist:
                    self.cxxclass.exclude(self, pid1, pid2)
                self.cxxclass.update(self)

    def exclude(self, pid1, pid2):
        if pmi.workerIsActive():
            self.cxxclass.exclude(self, pid1, pid2)

    def unexclude(self, pid1, pid2):
        if pmi.workerIsActive():
            self.cxxclass.unexclude(self, pid1, pid2)

    def observe_tuple(self, fpl):
        if pmi.workerIsActive():
            self.cxxclass.observe_tuple(self, fpl)

    def observe_triple(self, ftl):
        if pmi.workerIsActive():
            self.cxxclass.observe_triple(self, ftl)

    def observe_quadruple(self, fql):
        if pmi.workerIsActive():
            self.cxxclass.observe_quadruple(self, fql)

class VerletListLocal(_espressopp.VerletList):

    def __init__(self, system, cutoff, exclusionlist=None):
        if pmi.workerIsActive():
            if exclusionlist is None:
                # rebuild list in constructor
                cxxinit(self, _espressopp.VerletList, system, cutoff, True)
            elif isinstance(exclusionlist, DynamicExcludeListLocal):
                cxxinit(self, _espressopp.VerletList, system, cutoff, exclusionlist, True)
            else:
                # do not rebuild list in constructor
                cxxinit(self, _espressopp.VerletList, system, cutoff, False)
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
  class DynamicExcludeList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
        cls='espressopp.DynamicExcludeListLocal',
        pmiproperty=['is_dirty', 'size'],
        pmicall=['exclude', 'unexclude', 'connect', 'disconnect', 'observe_tuple',
                 'observe_triple', 'observe_quadruple', 'update'],
        pmiinvoke=['get_list']
    )

  class VerletList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls = 'espressopp.VerletListLocal',
      pmiproperty = [ 'builds' ],
      pmicall = [ 'totalSize', 'exclude', 'connect', 'disconnect', 'getVerletCutoff', 'setVerletCutoff' ],
      pmiinvoke = [ 'getAllPairs', 'get_timers', 'excludeListSize', 'rebuild']
    )
