#  Copyright (C) 2017
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
********************************************
**espressopp.integrator.LangevinThermostatOnRadius**
********************************************
Langevin Thermostat for Radii of Particles

Example:

>>> radius_mass = mass
>>> # set virtual mass for dynamics of radius
>>> langevin = espressopp.integrator.LangevinThermostatOnRadius(system, radius_mass)
>>> # set up the thermostat
>>> langevin.gamma = gamma
>>> # set friction coefficient gamma
>>> langevin.temperature = temp
>>> # set temperature
>>> integrator.addExtension(langevin)
>>> # add extensions to a previously defined integrator

.. function:: espressopp.integrator.LangevinThermostatOnRadius(system, dampingmass)

		:param system: 
		:param _dampingmass: 
		:type system: 
		:type dampingmass: real

.. function:: espressopp.integrator.LangevinThermostatOnRadius.addExclusions(pidlist)

                :param pidlist: list of particle ids to be excluded from thermostating.
                :type pidlist: list of ints
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_LangevinThermostatOnRadius 

class LangevinThermostatOnRadiusLocal(ExtensionLocal, integrator_LangevinThermostatOnRadius):

    def __init__(self, system, dampingmass):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LangevinThermostatOnRadius, system, dampingmass)

    def addExclusions(self, pidlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            for pid in pidlist:
                self.cxxclass.addExclpid(self, pid)

if pmi.isController :
    class LangevinThermostatOnRadius(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.LangevinThermostatOnRadiusLocal',
            pmiproperty = [ 'gamma', 'temperature'],
            pmicall = [ 'addExclusions' ]
            )
