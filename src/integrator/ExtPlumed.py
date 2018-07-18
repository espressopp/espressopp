#  Copyright (C) 2012,2013,2017
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
******************************
espressopp.integrator.ExtPlumed
******************************

usage:

plumed = espressopp.integrator.ExtPlumed(system, "plumed.dat", "log.plumed", 0.005)
plumed.setNaturalUnits()
plumed.Init()
integrator.addExtension(plumed)

or:

plumed = espressopp.integrator.ExtPlumed(system, "plumed.dat", "log.plumed", 0.005)
plumed.setNaturalUnits()
plumed.setRestart(1)
plumed.Init()
integrator.addExtension(plumed)

.. function:: espressopp.integrator.ExtPlumed(system, cmd, log, dt)

		:param system: The Espresso++ system object.
                :type system: espressopp.System
                :param cmd: input file for PLUMED
                :type cmd: ``str``
                :param log: log file for PLUMED
                :type log:  ``str``
                :param dt:  time step
                :type dt: ``float`` (default: 0.005)

"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_ExtPlumed
import mpi4py.MPI as MPI

class ExtPlumedLocal(ExtensionLocal, integrator_ExtPlumed):

    def __init__(self, system, cmd, log, dt):
        if pmi._PMIComm and pmi._PMIComm.isActive():
            if pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
                cxxinit(self, integrator_ExtPlumed, system, pmi._PMIComm.getMPIsubcomm(), cmd, log, dt)
            else:
                pass
        else:
            cxxinit(self, integrator_ExtPlumed, system, pmi._MPIcomm, cmd, log, dt)

if pmi.isController :
    class ExtPlumed(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.integrator.ExtPlumedLocal',
            pmicall = ['getBias', 'setNaturalUnits', 'setTimeUnit', 'setEnergyUnit', 'setLengthUnit', 'setKbT', 'setRealPrecison', 'setMDChargeUnit', 'setMDMassUnit', 'setRestart', 'readInputLine', 'Init', ]
            )
