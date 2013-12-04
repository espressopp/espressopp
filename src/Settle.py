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
*******************
**espresso.Settle**
*******************

"""
from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class SettleLocal(_espresso.Settle):
    'The (local) settle.'

    def __init__(self, storage, integrator, mO=16.0, mH=1.0, distHH=1.58, distOH=1.0):
        'Local construction of a settle class'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.Settle, storage, integrator, mO, mH, distHH, distOH)

    def addMolecules(self, moleculelist):
        """
        Each processor takes the broadcasted list.
        """
        if pmi.workerIsActive():
            for pid in moleculelist: 
                self.cxxclass.add(self, pid)


if pmi.isController:
    class Settle(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.SettleLocal',
            pmicall = [ "addMolecules" ]
            )
