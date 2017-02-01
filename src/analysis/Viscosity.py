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
*****************************
espressopp.analysis.Viscosity
*****************************


.. function:: espressopp.analysis.Viscosity(system)

		:param system:
		:type system:

.. function:: espressopp.analysis.Viscosity.compute(t0, dt, T)

		:param t0:
		:param dt:
		:param T:
		:type t0:
		:type dt:
		:type T:
		:rtype:

.. function:: espressopp.analysis.Viscosity.gather()

		:rtype:
"""
from espressopp import pmi
from espressopp.esutil import cxxinit

from _espressopp import analysis_Viscosity

from espressopp.analysis.Autocorrelation import *

class ViscosityLocal(AutocorrelationLocal, analysis_Viscosity):

    def __init__(self, system):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        cxxinit(self, analysis_Viscosity, system)

    def gather(self):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        return self.cxxclass.gather(self)

    def compute(self, t0, dt, T):
      if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
        return self.cxxclass.compute(self, t0, dt, T)

if pmi.isController:
  class Viscosity(Autocorrelation):

    #__metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.analysis.ViscosityLocal',
      pmicall = [ 'gather', 'compute' ]
    )
    def __init__(self, system):
      self.pmiinit(system)
