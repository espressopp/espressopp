#  Copyright (C) 2012,2013,2015
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
*******************************
espressopp.FixedTupleListAdress
*******************************

The FixedTupleListAdress is important for AdResS and H-AdResS simulations. It is the
connection between the atomistic and coarse-grained particles. It defines which
atomistic particles belong to which coarse-grained particle. In the following
example "tuples" is a python list of the form
( (pid_CG1, pidAT11,  pidAT12, pidAT13, ...), (pid_CG2, pidAT21,  pidAT22, pidAT23, ...), ...).
Each inner list (pid_CG1, pidAT11,  pidAT12, pidAT13, ...) defines a tuple. The
first number is the particle id of the coarse-grained particle while the
following numbers are the particle ids of the corresponding atomistic particles. 

Example - creating the FixedTupleListAdress:

>>> ftpl = espressopp.FixedTupleListAdress(system.storage)
>>> ftpl.addTuples(tuples)
>>> system.storage.setFixedTuples(ftpl)


.. function:: espressopp.FixedTupleListAdress(storage)

		:param storage: 
		:type storage: 

.. function:: espressopp.FixedTupleListAdress.addTuples(tuplelist)

		:param tuplelist: 
		:type tuplelist: 
		:rtype: 
"""

from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit

class FixedTupleListAdressLocal(_espressopp.FixedTupleListAdress):


    def __init__(self, storage):

        if pmi.workerIsActive():
            cxxinit(self, _espressopp.FixedTupleListAdress, storage)

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
                self.cxxclass.addTs(self);


if pmi.isController:
    class FixedTupleListAdress(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.FixedTupleListAdressLocal',
            pmicall = [ "addTuples" ]
            )
