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
**********************************************
espressopp.integrator.LangevinThermostatHybrid
**********************************************

As LangevinThermostat, but for use in AdResS systems, to allow the application of different thermostat friction constants (:math:`\gamma`) to different AdResS regions. Uses three values of :math:`\gamma`, one for the atomistic region, one for the hybrid region, and one for the coarse-grained region.

  >>> # create FixedTupleList object
  >>> ftpl = espressopp.FixedTupleListAdress(system.storage)
  >>> ftpl.addTuples(tuples)
  >>> system.storage.setFixedTuplesAdress(ftpl)
  >>>
  >>> system.storage.decompose()
  >>>
  >>> # create Langevin thermostat
  >>> thermostat             = espressopp.integrator.LangevinThermostatHybrid(system,ftpl)
  >>>
  >>> # set Langevin friction constants
  >>> thermostat.gamma       = 0.0 # units = 1/timeunit
  >>> print "# gamma for atomistic region for langevin thermostat = ",thermostat.gamma
  >>> thermostat.gammahy     = 10.0 # units = 1/timeunit
  >>> print "# gamma for hybrid region for langevin thermostat = ",thermostat.gammahy
  >>> thermostat.gammacg     = 10.0 # units = 1/timeunit
  >>> print "# gamma for coarse-grained region for langevin thermostat = ",thermostat.gammacg
  >>>
  >>> # set temperature of thermostat
  >>> thermostat.temperature = kBT
  >>> # kBT is a float with the value of temperature in reduced units, i.e. temperature * Boltzmann's constant in appropriate units

No need to include the line

  >>> thermostat.adress = True

as is necessary in the case of the basic LangevinThermostat, because LangevinThermostatHybrid is always only used in AdResS systems

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_LangevinThermostatHybrid

class LangevinThermostatHybridLocal(ExtensionLocal, integrator_LangevinThermostatHybrid):
    def __init__(self, system, fixedtuplelist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LangevinThermostatHybrid, system,fixedtuplelist)


if pmi.isController :
    class LangevinThermostatHybrid(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.LangevinThermostatHybridLocal',
            pmiproperty = [ 'gamma', 'gammahy','gammacg','temperature' ]
            )
