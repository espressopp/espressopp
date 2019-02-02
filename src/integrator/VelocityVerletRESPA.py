#  Copyright (C) 2018
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
***********************************************
espressopp.integrator.VelocityVerletRESPA
***********************************************

This is a multiple time stepping integrator according to the RESPA scheme (J. Chem. Phys. 97, 1990 (1992)). It has two layers: All forces of type "NonbondedSlow" are updated with a frequency given by the long time step, while all other forces are calculated according to the short time step. The short time step can be defined and set as a property of the integrator object, while the long time step is given by the product of the short time step with an integer "multistep", which can also be set.

Example:

>>> integrator = espressopp.integrator.VelocityVerletRESPA(system)
>>> integrator.dt = timestep
>>> integrator.multistep = multistep
>>> ...
>>> integrator.run(nsteps)

.. py:class:: espressopp.integrator.VelocityVerletRESPA(system)

        Constructs the VelocityVerletRESPA object.

        :param system: system object
        :type system: shared_ptr<System>

.. function:: espressopp.integrator.VelocityVerletRESPA.setmultistep(multistep)

        Sets the multiplier to construct the large timestep by multiplication with short time step as long_timestep = multistep * dt

        :param multistep: multiplier to construct the large timestep by multiplication with short time step
        :type sSteps: int

.. function:: espressopp.integrator.VelocityVerletRESPA.getmultistep()

        Gets the multiplier to construct the long timestep by multiplication with short time step as long_timestep = multistep * dt

        :return: multiplier to construct the long timestep by multiplication with short time step
        :rtype: int

.. py:data:: int espressopp.integrator.VelocityVerletRESPA.multistep

        Multiplier to construct the long timestep by multiplication with short time step as long_timestep = multistep * dt

.. py:data:: real espressopp.integrator.VelocityVerletRESPA.dt

        The short time step

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.MDIntegrator import *
from _espressopp import integrator_VelocityVerletRESPA

class VelocityVerletRESPALocal(MDIntegratorLocal, integrator_VelocityVerletRESPA):

    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_VelocityVerletRESPA, system)

    def setmultistep(self, multistep):
        if multistep <= 0:
            raise ValueError('multistep must be larger than zero. Your input: {}'.format(multistep))
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setmultistep(self, multistep)

    def getmultistep(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getmultistep(self)

if pmi.isController :
    class VelocityVerletRESPA(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
          cls =  'espressopp.integrator.VelocityVerletRESPALocal',
          pmiproperty = ['multistep'],
          pmicall = ['setmultistep', 'getmultistep']
        )
