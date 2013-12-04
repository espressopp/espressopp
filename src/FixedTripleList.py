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
****************************
**espresso.FixedTripleList**
****************************

"""
from espresso import pmi
import _espresso
import espresso
from espresso.esutil import cxxinit

class FixedTripleListLocal(_espresso.FixedTripleList):
    'The (local) fixed triple list.'

    def __init__(self, storage):
        'Local construction of a fixed triple list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedTripleList, storage)

    def add(self, pid1, pid2, pid3):
        'add triple to fixed triple list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3)

    def addTriples(self, triplelist):
        """
        Each processor takes the broadcasted triplelist and
        adds those triples whose first particle is owned by
        this processor.
        """
        if pmi.workerIsActive():
            for triple in triplelist:
                pid1, pid2, pid3 = triple
                self.cxxclass.add(self, pid1, pid2, pid3)

    def size(self):
        'count number of Triples in GlobalTripleList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    '''
    def addTriples(self, triplelist):
        """
        Each processor takes the broadcasted triplelist and
        adds those triples whose first particle is owned by
        this processor.
        """
        
        if pmi.workerIsActive():
            for triple in triplelist:
                pid1, pid2, pid3 = triple
                self.cxxclass.add(self, pid1, pid2, pid3)
    '''

    def getTriples(self):
        'return the triples of the GlobalTripleList'
        if pmi.workerIsActive():
          triples = self.cxxclass.getTriples(self)
          return triples 

if pmi.isController:
    class FixedTripleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedTripleListLocal',
            localcall = [],
            pmicall = [ "add", "addTriples" ],
            pmiinvoke = ["getTriples", "size"]
        )
