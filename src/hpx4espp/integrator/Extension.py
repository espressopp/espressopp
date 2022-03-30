#  Copyright (C) 2020-2022
#      Max Planck Institute for Polymer Research & JGU Mainz
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
****************************************
espressopp.hpx4espp.integrator.Extension
****************************************



.. function:: espressopp.hpx4espp.integrator.Extension.connect()

		:rtype:

.. function:: espressopp.hpx4espp.integrator.Extension.disconnect()

		:rtype:
"""
#from espressopp.esutil import cxxinit
from espressopp import pmi
from _espressopp import integrator_Extension

class ExtensionLocal(object):

    def connect(self):
        return self.cxxclass.connect(self)

    def disconnect(self):
        return self.cxxclass.disconnect(self)

if pmi.isController :
    class Extension(object):
        pmiproxydefs = dict(
            #cls =  'espressopp.hpx4espp.integrator.Extension',
            pmiproperty = [ 'type'],
            #pmicall = ['addForce']
            pmicall = [ 'connect', 'disconnect' ]
        )
