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
***************************************
espressopp.integrator.BerendsenBarostat
***************************************

This is the Berendsen barostat implementation according to the original paper [Berendsen84]_.
If Berendsen barostat is defined (as a property of integrator) then at the each run the system size
and the particle coordinates will be scaled by scaling parameter :math:`\mu` according to
the formula:

.. math::   \mu = [1 - \Delta t/\tau (P_{0} - P)]^{1/3}

where :math:`\Delta t` - integration timestep, :math:`\tau` - time parameter (coupling parameter),
:math:`P_{0}` - external pressure and :math:`P` - instantaneous pressure. 

Example:

    >>> berendsenP = espressopp.integrator.BerendsenBarostat(system)
    >>> berendsenP.tau = 0.1
    >>> berendsenP.pressure = 1.0
    >>> integrator.addExtension(berendsenP)

**!IMPORTANT** In order to run *npt* simulation one should separately define thermostat as well 
(e.g. BerendsenThermostat_).

.. _BerendsenThermostat: espressopp.integrator.BerendsenThermostat.html

Definition:

    In order to define the Berendsen barostat
    
    >>> berendsenP = espressopp.integrator.BerendsenBarostat(system)
    
    one should have the System_ defined.

.. _System: espressopp.System.html

Properties:

*   *berendsenP.tau*

    The property 'tau' defines the time parameter :math:`\tau`.

*   *berendsenP.pressure*
    
    The property 'pressure' defines the external pressure :math:`P_{0}`.
    
Setting the integration property:
    
    >>> integrator.addExtension(berendsenP)
    
    It will define Berendsen barostat as a property of integrator.
    
One more example:

    >>> berendsen_barostat = espressopp.integrator.BerendsenBarostat(system)
    >>> berendsen_barostat.tau = 10.0
    >>> berendsen_barostat.pressure = 3.5
    >>> integrator.addExtension(berendsen_barostat)


Canceling the barostat:
    
    If one do not need the pressure regulation in system anymore or need to switch
    the ensamble or whatever :)

    >>> # define barostat with parameters
    >>> berendsen = espressopp.integrator.BerendsenBarostat(system)
    >>> berendsen.tau = 0.8
    >>> berendsen.pressure = 15.0
    >>> integrator.addExtension(berendsen)
    >>> ...
    >>> # some runs
    >>> ...
    >>> # disconnect Berendsen barostat
    >>> berendsen.disconnect()
    >>> # the next runs will not include the system size and particle coordinates scaling

    
  Connecting the barostat back after the disconnection
  
    >>> berendsen.connect()

References:

.. [Berendsen84] Berendsen et al., *J. Chem. Phys.*, 81, **1984**, p. 3684


.. function:: espressopp.integrator.BerendsenBarostat(system)

		:param system: 
		:type system: 
"""


from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_BerendsenBarostat

class BerendsenBarostatLocal(ExtensionLocal, integrator_BerendsenBarostat):
  def __init__(self, system):

    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or \
            pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, integrator_BerendsenBarostat, system)

if pmi.isController:
  class BerendsenBarostat(Extension):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.integrator.BerendsenBarostatLocal',
      pmiproperty = [ 'tau', 'pressure', 'fixed' ]
    )
