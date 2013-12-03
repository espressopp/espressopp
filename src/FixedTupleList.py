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
***************************
**espresso.FixedTupleList**
***************************

"""
from espresso import pmi
import _espresso
import espresso
from espresso.esutil import cxxinit

class FixedTupleListLocal(_espresso.FixedTupleList):
    'The (local) fixed tuple list.'

    def __init__(self, storage):
        'Local construction of a fixed tuple list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedTupleList, storage)

    """def addTuples(self, tuplelist):
        'add tuple to fixed tuple list'
        if pmi.workerIsActive():
            return self.cxxclass.addTuple(self, tuplelist)"""


    def size(self):
        'count number of Tuple in GlobalTupleList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.size(self)



if pmi.isController:
    class FixedTupleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedTupleListLocal',
            #localcall = [ "add" ],
            pmicall = [ "addTuple", "getTuples" ],
            pmiinvoke = ["size"]
        )
