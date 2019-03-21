#  Copyright (C) 2018, 2017
#      Max Planck Institute for Polymer Research
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
**********************************
**espressopp.FixedLocalTupleList**
**********************************
This class can contain many tuple which store a arbitrary positive number, which should be more than 2, of local (real + ghost) particle id.

For using this class, there is 1 conditions:

Particles in one tuple must be in a same or neighbor cell list.

.. function:: espressopp.FixedLocalTupleList(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.FixedLocalTupleList.addTuple(tuple)

		:param tuple: 
		:type tuple: python::list 


.. function:: espressopp.FixedLocalTupleList.getTuples()

		:rtype: 

.. function:: espressopp.FixedLocalTupleList.size()

		:rtype: 
"""
from espressopp import pmi
import _espressopp
import espressopp
from espressopp.esutil import cxxinit

class FixedLocalTupleListLocal(_espressopp.FixedLocalTupleList):


    def __init__(self, storage):

        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedLocalTupleList, storage)

    def getTuples(self):

        if pmi.workerIsActive():
            tuples=self.cxxclass.getTuples(self)
            return tuples

    def size(self):

        if pmi.workerIsActive():
            return self.cxxclass.size(self)



if pmi.isController:
    class FixedLocalTupleList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.FixedLocalTupleListLocal',
            #localcall = [ "add" ],
            pmicall = [ "addTuple"],
            pmiinvoke = ["size", "getTuples"]
        )
