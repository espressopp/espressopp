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
***************************
espressopp.VerletListTriple
***************************


.. function:: espressopp.VerletListTriple(system, cutoff, exclusionlist)

		:param system: 
		:param cutoff: 
		:param exclusionlist: (default: [])
		:type system: 
		:type cutoff: 
		:type exclusionlist: 

.. function:: espressopp.VerletListTriple.exclude(exclusionlist)

		:param exclusionlist: 
		:type exclusionlist: 
		:rtype: 

.. function:: espressopp.VerletListTriple.getAllTriples()

		:rtype: 

.. function:: espressopp.VerletListTriple.localSize()

		:rtype: 

.. function:: espressopp.VerletListTriple.totalSize()

		:rtype: 
"""
from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit

class VerletListTripleLocal(_espressopp.VerletListTriple):


    def __init__(self, system, cutoff, exclusionlist=[]):

        if pmi.workerIsActive():
          '''
          cxxinit(self, _espressopp.VerletListTriple, system, cutoff, True)
          if (exclusionlist != []):
            print 'Warning! Exclusion list is not yet implemented to the triple verlet \
                  list. Nothing happend to exclusion list'
          '''

          if (exclusionlist == []):
            # rebuild list in constructor
            cxxinit(self, _espressopp.VerletListTriple, system, cutoff, True)
          else:
            # do not rebuild list in constructor
            cxxinit(self, _espressopp.VerletListTriple, system, cutoff, False)
            # add exclusions
            for pid in exclusionlist:
                self.cxxclass.exclude(self, pid)
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
        #if pmi.workerIsActive():
        #  print 'Warning! Exclusion list is not yet implemented to the triple verlet \
        #        list. Nothing happend'
        for pid in exclusionlist:
            self.cxxclass.exclude(self, pid)
        # rebuild list with exclusions
        self.cxxclass.rebuild(self)
            
    def getAllTriples(self):

        if pmi.workerIsActive():
            triples=[]
            ntriples=self.localSize()
            for i in xrange(ntriples):
              triple=self.cxxclass.getTriple(self, i+1)
              triples.append(triple)
            return triples


if pmi.isController:
  class VerletListTriple(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls = 'espressopp.VerletListTripleLocal',
      pmiproperty = [ 'builds' ],
      pmicall = [ 'totalSize', 'exclude', 'connect', 'disconnect', 'getVerletCutoff' ],
      pmiinvoke = [ 'getAllTriples' ]
    )
