#  Copyright (C) 2012,2013,2014,2015,2016
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
espressopp.integrator.LangevinThermostat
****************************************

Langevin Thermostat

Example:

>>> langevin = espressopp.integrator.LangevinThermostat(system)
>>> # set up the thermostat
>>> langevin.gamma = gamma
>>> # set friction coefficient gamma
>>> langevin.temperature = temp
>>> # set temperature
>>> langevin.adress = True
>>> # set adress (default is False)
>>> integrator.addExtension(langevin)
>>> # add extensions to a previously defined integrator

.. function:: espressopp.integrator.LangevinThermostat(system)

        :param system: system object
        :type system: shared_ptr<System>

.. function:: espressopp.integrator.LangevinThermostat.addExclusions(pidlist)

        :param pidlist: list of particle ids to be excluded from thermostating. In adaptive (AdResS) simulations, add ids of atomistic particles to be excluded (thermostats acts in this case on atomistic level). For normal simulations, add normal or coarse-grained particle ids.
        :type pidlist: list of ints

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_LangevinThermostat

class LangevinThermostatLocal(ExtensionLocal, integrator_LangevinThermostat):

    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_LangevinThermostat, system)

    def addExclusions(self, pidlist):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            for pid in pidlist:
                self.cxxclass.addExclpid(self, pid)

if pmi.isController :
    class LangevinThermostat(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.LangevinThermostatLocal',
            pmiproperty = [ 'gamma', 'temperature', 'adress' ],
            pmicall = [ 'addExclusions' ]
            )
