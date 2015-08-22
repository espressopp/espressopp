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


r"""
**************************************
**AdResS** - Object
**************************************

The AdResS object is an extension to the integrator. It makes sure that the
integrator also processes the atomistic particles and not only the CG particles.
Hence, this object is of course only used when performing AdResS or H-AdResS
simulations.

In detail the AdResS extension makes sure:
---------------------------------------------

* that also the forces on the atomistic particles are initialized and set to
  by Adress::initForces
* that also the atomistic particles are integrated and propagated by
  Adress::integrate1 and Adress::integrate2

Example - how to turn on the AdResS integrator extension:

>>> adress      = espressopp.integrator.Adress(system)
>>> integrator.addExtension(adress)


.. function:: espressopp.integrator.Adress(_system, _verletlist, _fixedtuplelist, KTI)

		:param _system: 
		:param _verletlist: 
		:param _fixedtuplelist: 
		:param KTI: (default: False)
		:type _system: 
		:type _verletlist: 
		:type _fixedtuplelist: 
		:type KTI: 
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_Adress

class AdressLocal(ExtensionLocal, integrator_Adress):


    def __init__(self, _system, _verletlist, _fixedtuplelist, KTI = False):

        if pmi.workerIsActive():
            cxxinit(self, integrator_Adress, _system, _verletlist, _fixedtuplelist, KTI)

if pmi.isController:
    class Adress(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.integrator.AdressLocal' #,
            #pmiproperty = [ 'builds' ],
            #pmicall = [ 'totalSize', 'exclude', 'addAdrParticles', 'rebuild' ]
            )
