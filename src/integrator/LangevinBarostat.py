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
**************************************
espressopp.integrator.LangevinBarostat
**************************************

This is the barostat implementation to perform Langevin dynamics in a Hoover style extended system
according to the paper [Quigley04]_. It includes corrections of Hoover approach which were introduced
by Martyna et al [Martyna94]_.
If LangevinBarostat is defined (as a property of integrator) the integration equations will be 
modified. The volume of system :math:`V` is introduced as a dynamical variable:

.. math:: \boldsymbol{\dot{r}}_{i} = \frac{\boldsymbol{p}_{i}}{m_{i}} +
          \frac{p_{\epsilon}}{W}\boldsymbol{r}_{i}
					
.. math:: \boldsymbol{\dot{p}}_{i} = -\bigtriangledown_{\boldsymbol{r}_{i}}\Phi -
          (1+\frac{n}{N_{f}})\frac{p_{\epsilon}}{W}\boldsymbol{p}_{i} - 
          \gamma\boldsymbol{p}_{i} + \boldsymbol{R}_{i}
					
.. math:: \dot{V} = dVp_{\epsilon}/W

.. math:: \dot{p}_{\epsilon} = nV(X-P_{ext})+
          \frac{n}{N_{f}}\sum^{N}_{i=1}\frac{\boldsymbol{p}_{i}^{2}}{m_{i}} -
          \gamma_{p}p_{\epsilon} + R_{p}


where volume has a fictitious mass :math:`W` and associated momentum :math:`p_{\epsilon}`,
:math:`\gamma_{p}` - friction coefficient,
:math:`P_{ext}` - external pressure and :math:`X` - instantaneous pressure without white noise
contribution from thermostat, :math:`n` - dimension, :math:`N_{f}` - degrees of freedom (if there
are no constrains and :math:`N` is the number of particles in system :math:`N_{f}=nN`).
:math:`R_{p}` - values which are drawn from Gaussian distribution of zero mean and unit variance
scaled by

.. math:: \sqrt{\frac{2k_{B}TW\gamma_{p}}{\Delta t}}

**!IMPORTANT** Terms :math:`- \gamma\boldsymbol{p}_{i} + \boldsymbol{R}_{i}` correspond to the 
termostat. They are not included here and will not be calculated if the Langevin Thermostat is not
defined.

Example:

    >>> rng = espressopp.esutil.RNG()
    >>> langevinP = espressopp.integrator.LangevinBarostat(system, rng, desiredTemperature)
    >>> langevinP.gammaP = 0.05
    >>> langevinP.pressure = 1.0
    >>> langevinP.mass = pow(10.0, 4)
    >>> integrator.addExtension(langevinP)

**!IMPORTANT**  This barostat is supposed to be run in a couple with thermostat in order 
to simulate the *npt* ensamble, because the term :math:`R_{p}` needs the temperature as a 
parameter. 

Definition:

    In order to define the Langevin-Hoover barostat
    
    >>> langevinP = espressopp.integrator.LangevinBarostat(system, rng, desiredTemperature)
    
    one should have the System_ and RNG_ defined and know the desired temperature.

.. _System: espressopp.System.html
.. _RNG:

Properties:

*   *langevinP.gammaP*

    The property 'gammaP' defines the friction coefficient :math:`\gamma_{p}`.

*   *langevinP.pressure*
    
    The property 'pressure' defines the external pressure :math:`P_{ext}`.
    
*   *langevinP.mass*

    The property 'mass' defines the fictitious mass :math:`W`.

Methods:

*   *setMassByFrequency( frequency )*

    Set the proper *langevinP.mass* using expression :math:`W=dNk_{b}T/\omega_{b}^{2}`, where
    frequency, :math:`\omega_{b}`, is the frequency of required volume fluctuations. The value of
    :math:`\omega_{b}` should be less then the lowest frequency which appears in the NVT temperature
    spectrum [Quigley04]_ in order to match the canonical distribution. :math:`d` - dimensions,
    :math:`N` - number of particles, :math:`k_{b}` - Boltzmann constant, :math:`T` - desired
    temperature.
    
**NOTE** The *langevinP.mass* can be set both directly and using the 
(*setMassByFrequency( frequency )*) 

Adding to the integration:
    
    >>> integrator.addExtension(langevinP)
    
    It will define Langevin-Hoover barostat as a property of integrator.
    
One more example:

    >>> rngBaro = espressopp.esutil.RNG()
    >>> lP = espressopp.integrator.LangevinBarostat(system, rngBaro, desiredTemperature)
    >>> lP.gammaP = .5
    >>> lP.pressure = 1.0
    >>> lP.mass = pow(10.0, 5)
    >>> integrator.addExtension(lP)

Canceling the barostat:
    
    If one do not need the pressure regulation in system anymore or need to switch
    the ensamble or whatever :)

    >>> # define barostat with parameters
    >>> rngBaro = espressopp.esutil.RNG()
    >>> lP = espressopp.integrator.LangevinBarostat(system, rngBaro, desiredTemperature)
    >>> lP.gammaP = .5
    >>> lP.pressure = 1.0
    >>> lP.mass = pow(10.0, 5)
    >>> integrator.langevinBarostat = lP
    >>> ...
    >>> # some runs
    >>> ...
    >>> # disconnect barostat
    >>> langevinBarostat.disconnect()
    >>> # the next runs will not include the modification of integration equations

    
  Connecting the barostat back after the disconnection
  
    >>> langevinBarostat.connect()

References:

.. [Quigley04] D. Quigley, M.I.J. Probert, *J. Chem. Phys.*, 120, **2004**, p. 11432

.. [Martyna94] G. Martyna, D. Tobias, M. Klein, *J. Chem. Phys.*, 101, **1994**, p. 4177


.. function:: espressopp.integrator.LangevinBarostat(system, rng, temperature)

		:param system: 
		:param rng: 
		:param temperature: 
		:type system: 
		:type rng: 
		:type temperature: 
"""


from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_LangevinBarostat

class LangevinBarostatLocal(ExtensionLocal, integrator_LangevinBarostat):
  def __init__(self, system, rng, temperature):

    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, integrator_LangevinBarostat, system, rng, temperature)

if pmi.isController :
  class LangevinBarostat(Extension):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.integrator.LangevinBarostatLocal',
      pmiproperty = [ 'gammaP', 'pressure', 'mass' ],
      pmicall = [ "setMassByFrequency" ]
    )
