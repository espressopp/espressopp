#  Copyright (C) 2016,2021
#      Jakub Krajniak (jkrajniak at gmail.com)
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
**************************************
**FixedVSList** - Object
**************************************

The FixedVSList keeps the connection between the coarse-grained and underlaying atomistic
particles.

Example - creating the FixedVSList:

>>> tuples = [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
>>> ftpl = espressopp.FixedVSList(system.storage)
>>> ftpl.addTuples(tuples)



.. function:: espressopp.FixedVSList(storage)
		:param storage:
		:type storage: 

.. function:: espressopp.FixedVSList.addTuples(tuplelist)

		:param tuplelist: 
		:type tuplelist: 
		:rtype: 
"""

from espressopp import pmi
import _espressopp
import espressopp
from espressopp.esutil import cxxinit


class FixedVSListLocal(_espressopp.FixedVSList):
    def __init__(self, storage):
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedVSList, storage)

    def addTuples(self, tuplelist):
        """
        Each processor takes the broadcasted tuplelist and
        adds those tuples whose virtual particle is owned by
        this processor.
        """
        if pmi.workerIsActive():
            for tuple in tuplelist:
                for pid in tuple:
                    self.cxxclass.add(self, pid)
                self.cxxclass.addTs(self)


if pmi.isController:
    class FixedVSList(metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls='espressopp.FixedVSListLocal',
            pmicall=["addTuples"]
        )
