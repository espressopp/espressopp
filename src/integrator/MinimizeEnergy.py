#  Copyright (C) 2016
#      Jakub Krajniak (jkrajniak at gmail.com)
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
**espressopp.integrator.MinimizeEnergy**
****************************************

This is a very simple approach to perform energy minimialization of the system. The module uses
a `steepest descent method <https://en.wikipedia.org/wiki/Gradient_descent>`_. The position of particles is
updated following the equation:

.. math::

   p_{i+1} = p_i + min(\gamma F_i, d_{max})

where :math:`p_{i+}` is a new position, :math:`p_i` is a position at current step with corresponding
force :math:`F_i`. The parameters :math:`\gamma` and :math:`d_{max}` are set by user and control the relaxation
of the energy and the maximum update of the coordinates per step.
The routine runs until the maximal force is bigger than :math:`f_{max}` or for at most *n* steps.

**Please note**
This module does not support any integrator extensions.

.. function:: espressopp.integrator.MinimizeEnergy(system, gamma, max_force, max_displacement)

		:param system: The espressopp system object.
		:type system: espressopp.System
		:param gamma: The gamma value.
		:type gamma: float
		:param max_force: The maximum force threshold.
		:type max_force: float
		:param max_displacement: The maximum displacemenet.
		:type max_displacement: float

.. function:: espressopp.integrator.MinimizeEnergy.run(max_steps)

        :param max_steps: The maximum number of steps to run.
        :type max_steps: int

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from _espressopp import integrator_MinimizeEnergy

class MinimizeEnergyLocal(integrator_MinimizeEnergy):
    def __init__(self, system, gamma, max_force, max_displacement):
        if pmi.workerIsActive():
            cxxinit(self, integrator_MinimizeEnergy, system, gamma, max_force, max_displacement)

    def run(self, niter):
        if pmi.workerIsActive():
            return self.cxxclass.run(self, niter)

if pmi.isController:
    class MinimizeEnergy:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.MinimizeEnergyLocal',
            pmiproperty = ('f_max', 'displacement', 'step'),
            pmicall = ('run', )
        )
