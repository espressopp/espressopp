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
**********************************
**espresso.FixedTripleListAdress**
**********************************

"""
from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class FixedTripleListAdressLocal(_espresso.FixedTripleListAdress):
    'The (local) fixed triple list.'

    def __init__(self, storage, fixedtupleList):
        'Local construction of a fixed pair list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedTripleListAdress, storage, fixedtupleList)

    def add(self, pid1, pid2):
        'add pair to fixed triple list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2, pid3)

    def addTriples(self, triplelist):
        """
        Each processor takes the broadcasted triplelist and
        adds those pairs whose first particle is owned by
        this processor.
        """
        
        if pmi.workerIsActive():
            for triple in triplelist:
                pid1, pid2, pid3 = triple
                self.cxxclass.add(self, pid1, pid2, pid3)

if pmi.isController:
    class FixedTripleListAdress(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedTripleListAdressLocal',
            localcall = [ "add" ],
            pmicall = [ "addTriples" ]
            )
