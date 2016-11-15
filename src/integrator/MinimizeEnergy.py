#  Copyright (C) 2016
#      Max Planck Institute for Polymer Research
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
************************************
espressopp.integrator.MinimizeEnergy
************************************

This is a very simple approach to perform energy minimization of the system. The module uses
a `steepest descent method <https://en.wikipedia.org/wiki/Gradient_descent>`_. The position of particles is updated following the equation:

.. math::

   p_{i+1} = p_i + min(\gamma F_i, d_{max})

where :math:`p_{i+}` is a new position, :math:`p_i` is a position at current step with corresponding
force :math:`F_i`. The parameters :math:`\gamma` and :math:`d_{max}` are set by user and control the relaxation
of the energy and the maximum update of the coordinates per step.

Additionaly, a variable :math:`\gamma` step is also implemented. In this case, the position of particles is updated following the equation:

.. math::

   p_{i+1} = p_i + d_{max}/f_{max} F_i

where :math:`f_{max}` is a maximum force in a single step of steepest descent method. :math:`\gamma=d_{max}/f_{max}` is automatically adjusted to a force magnitude.

In both cases, the routine runs until the maximum force is bigger than :math:`f_{max}` or for at most *n* steps.

**Please note**
This module does not support any integrator extensions.

Example

>>> em = espressopp.integrator.MinimizeEnergy(system, gamma=0.001, ftol=0.01, max_displacement=0.0001)
>>> em.run(10000)

Example

>>> em = espressopp.integrator.MinimizeEnergy(system, gamma=0.01, ftol=0.01, max_displacement=0.01, variable_step_flag=True)
>>> em.run(10000)

**API**

.. function:: espressopp.integrator.MinimizeEnergy(system, gamma, ftol, max_displacement, variable_step_flag)

		:param system: The espressopp system object.
		:type system: espressopp.System
		:param gamma: The gamma value.
		:type gamma: float
		:param ftol: The force tolerance
		:type ftol: float
		:param max_displacement: The maximum displacement.
		:type max_displacement: float
                :param variable_step_flag: The flag of adjusting gamma to the force strength.
		:type variable_step_flag: bool

.. function:: espressopp.integrator.MinimizeEnergy.run(max_steps, verbose)

        :param max_steps: The maximum number of steps to run.
        :type max_steps: int
        :param verbose: If set to True then display information about maximum force during the iterations.
        :type verbose: bool
        :return: The true if the maximum force in the system is lower than ftol otherwise false.
        :rtype: bool

.. py:data:: f_max

    The maximum force in the system.

.. py:data:: displacement

    The maximum displacement used during the run of MinimizeEnergy

.. py:data:: step

    The current iteration step.

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from _espressopp import integrator_MinimizeEnergy

class MinimizeEnergyLocal(integrator_MinimizeEnergy):
    def __init__(self, system, gamma, ftol, max_displacement, variable_step_flag=False):
        if pmi.workerIsActive():
            cxxinit(self, integrator_MinimizeEnergy, system, gamma, ftol*ftol, max_displacement, variable_step_flag)

    def run(self, niter, verbose=False):
        if pmi.workerIsActive():
            return self.cxxclass.run(self, niter, verbose)

if pmi.isController:
    class MinimizeEnergy:
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.MinimizeEnergyLocal',
            pmiproperty = ('f_max', 'displacement', 'step'),
            pmicall = ('run', )
        )
