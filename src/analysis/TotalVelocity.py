#  Copyright (C) 2014 Pierre de Buyl
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
*************************************
**espressopp.analysis.TotalVelocity**
*************************************


.. function:: espressopp.analysis.TotalVelocity(system)

		:param system: The system object.
		:type system: espressopp.System

.. function:: espressopp.analysis.TotalVelocity.compute()

        Compute the total velocity of the system.

		:rtype: float

.. function:: espressopp.analysis.TotalVelocity.reset()

        Subtract the total velocity of the system from every particle.

Examples
---------

Reset the velocity
+++++++++++++++++

>>> total_velocity = espressopp.analysis.TotalVelocity(system)
>>> total_velocity.reset()

Extension to integrator
++++++++++++++++++++++++++++++++++++++++++++

This extension can also be attached to integrator and run `reset()` every `n-th` steps.

>>> total_velocity = espressopp.analysis.TotalVelocity(system)
>>> ext_remove_com = espressopp.analysis.ExtAnalyze(total_velocity, 10)
>>> integrator.addExtension(ext_remove_com)

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_TotalVelocity

class TotalVelocityLocal(ObservableLocal, analysis_TotalVelocity):

    def __init__(self, system):
	if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
          cxxinit(self, analysis_TotalVelocity, system)
    def compute(self):
        return self.cxxclass.compute(self)
    def reset(self):
        return self.cxxclass.reset(self)

if pmi.isController :
    class TotalVelocity(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.TotalVelocityLocal',
            pmicall = [ "compute", "reset" ],
            pmiproperty = ["v"]
            )
