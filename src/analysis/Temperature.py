#  Copyright (C) 2017
#      Jakub Krajniak (jkrajniak at gmail.com)
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
*******************************
espressopp.analysis.Temperature
*******************************

Calculate the temperature of the system (in :math:`k_B T` units).

.. function:: espressopp.analysis.Temperature(system)

    :param shared_ptr system: system object
    :returns: temperature
    :rtype: real

    Temperature of the system of :math:`N` particles is calculated as:

    .. math::

        T = \frac{1}{N_f} \sum^N_{i=1} m_i v_i^2,

    where :math:`m_i` and :math:`v_i` are the mass and velocity of a
    particle :math:`i`.

    :math:`N_f = 3N` is the number of the system's degrees of freedom.

    **Example:**

    >>> # declare an object, e.g., T:
    >>> T   = espressopp.analysis.Temperature(system)
    >>>
    >>> # later in your script compute temperature and print it:
    >>> print T.compute()


.. function::espressopp.analysis.Temperature.add_type(type_id)

      :param type_id: The particle type id to observe

.. function::espressopp.analysis.Temperature.remove_type(type_id)

      :param type_id: Remove particle type id.

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_Temperature

class TemperatureLocal(ObservableLocal, analysis_Temperature):

    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_Temperature, system)
    
    def add_type(self, type_id):
        if pmi.workerIsActive():
            self.cxxclass.add_type(self, type_id)

    def remove_type(self, type_id):
        if pmi.workerIsActive():
            ret_val = self.cxxclass.remove_type(self, type_id)
            if not ret_val:
                print('Warning, type {} not found'.format(type_id))


if pmi.isController :
    class Temperature(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.TemperatureLocal',
            pmicall = ['add_type', 'remove_type']
            )
