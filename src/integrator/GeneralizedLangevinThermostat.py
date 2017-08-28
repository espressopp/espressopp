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
***************************************************
espressopp.integrator.GeneralizedLangevinThermostat
***************************************************


.. function:: espressopp.integrator.GeneralizedLangevinThermostat(system)

		:param system: 
		:type system: 

.. function:: espressopp.integrator.GeneralizedLangevinThermostat.addCoeffs(itype, filename, type)

		:param itype: 
		:param filename: 
		:param type: 
		:type itype: 
		:type filename: 
		:type type: 
		:rtype: 
"""

# NOTE: This Generalized Langevin Thermostat uses the extended variable scheme
# and only applies the memory-based friction term, no noise is applied.

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_GeneralizedLangevinThermostat 

class GeneralizedLangevinThermostatLocal(ExtensionLocal, integrator_GeneralizedLangevinThermostat):

    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_GeneralizedLangevinThermostat, system)


    def addCoeffs(self, itype, filename, type):
            """
            Each processor takes the broadcasted interpolation type,
            filename and particle type
            """
            if pmi.workerIsActive():
                self.cxxclass.addCoeffs(self, itype, filename, type)
    #def enableAdress(self):
    #    if pmi.workerIsActive():
    #        self.cxxclass.enableAdress(self);

if pmi.isController :
    class GeneralizedLangevinThermostat(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.GeneralizedLangevinThermostatLocal',
            pmiproperty = [ 'itype', 'filename'],
            pmicall = ['addCoeffs' ]
            )
