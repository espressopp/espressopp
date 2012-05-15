"""
******************************************************
**LangevinBarostat** - Langevin-Hoover barostat Object
******************************************************

This is the barostat implementation to perform Langevin dynamics in a Hoover style extended system
according to the paper [Quigley04]_. It includes corrections of Hoover approach which were introduced
by Martyna et al [Martyna94]_.
If LangevinBarostat is defined (as a property of integrator) the integration equations will be 
modified. The volume of system :math:`V` is introduced as a dynamical variable:

.. math:: \\boldsymbol{\dot{r}}_{i} = \\frac{\\boldsymbol{p}_{i}}{m_{i}} +
          \\frac{p_{\epsilon}}{W}\\boldsymbol{r}_{i}
.. math:: \\boldsymbol{\dot{p}}_{i} = -\\bigtriangledown_{\\boldsymbol{r}_{i}}\Phi -
          (1+\\frac{n}{N_{f}})\\frac{p_{\epsilon}}{W}\\boldsymbol{p}_{i} - 
          \gamma\\boldsymbol{p}_{i} + \\boldsymbol{R}_{i}
.. math:: \\dot{V} = dVp_{\epsilon}/W
.. math:: \\dot{p}_{\epsilon} = nV(X-P_{ext})+
          \\frac{n}{N_{f}}\sum^{N}_{i=1}\\frac{\\boldsymbol{p}_{i}^{2}}{m_{i}} -
          \gamma_{p}p_{\epsilon} + R_{p}


where volume has a fictitious mass :math:`W` and associated momentum :math:`p_{\epsilon}`,
:math:`\gamma_{p}` - friction coefficient,
:math:`P_{ext}` - external pressure and :math:`X` - instantaneous pressure without white noise
contribution from thermostat, :math:`n` - dimension, :math:`N_{f}` - degrees of freedom (if there
are no constrains and :math:`N` is the number of particles in system :math:`N_{f}=nN`).
:math:`R_{p}` - values which are drawn from Gaussian distribution of zero mean and unit variance
scaled by

.. math:: \sqrt{\\frac{2k_{B}TW\gamma_{p}}{\Delta t}}

**!IMPORTANT** Terms :math:`- \gamma\\boldsymbol{p}_{i} + \\boldsymbol{R}_{i}` correspond to the 
termostat. They are not included here and will not be calculated if the Langevin Thermostat is not
defined.

Example:

    >>> rng = espresso.esutil.RNG()
    >>> langevinP = espresso.integrator.LangevinBarostat(system, rng, desiredTemperature)
    >>> langevinP.gammaP = 0.05
    >>> langevinP.pressure = 1.0
    >>> langevinP.mass = pow(10.0, 4)
    >>> integrator.langevinBarostat = langevinP

**!IMPORTANT**  This barostat is supposed to be run in a couple with thermostat in order 
to simulate the *npt* ensamble, because the term :math:`R_{p}` needs the temperature as a 
parameter. 

Definition:

    In order to define the Langevin-Hoover barostat
    
    >>> langevinP = espresso.integrator.LangevinBarostat(system, rng, desiredTemperature)
    
    one should have the System_ and RNG_ defined and know the desired temperature.

.. _System: espresso.System.html
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

Setting the integration property:
    
    >>> integrator.langevinBarostat = langevinP
    
    It will define Langevin-Hoover barostat as a property of integrator.
    
One more example:

    >>> rngBaro = espresso.esutil.RNG()
    >>> lP = espresso.integrator.LangevinBarostat(system, rngBaro, desiredTemperature)
    >>> lP.gammaP = .5
    >>> lP.pressure = 1.0
    >>> lP.mass = pow(10.0, 5)
    >>> integrator.langevinBarostat = lP

Canceling the barostat:
    
    If one do not need the pressure regulation in system anymore or need to switch
    the ensamble or whatever :)

    >>> # define barostat with parameters
    >>> rngBaro = espresso.esutil.RNG()
    >>> lP = espresso.integrator.LangevinBarostat(system, rngBaro, desiredTemperature)
    >>> lP.gammaP = .5
    >>> lP.pressure = 1.0
    >>> lP.mass = pow(10.0, 5)
    >>> integrator.langevinBarostat = lP
    >>> ...
    >>> # some runs
    >>> ...
    >>> # erase barostat
    >>> integrator.langevinBarostat = None
    >>> # the next runs will not include the modification of integration equations

References:

.. [Quigley04] D. Quigley, M.I.J. Probert, *J. Chem. Phys.*, 120, **2004**, p. 11432

.. [Martyna94] G. Martyna, D. Tobias, M. Klein, *J. Chem. Phys.*, 101, **1994**, p. 4177

"""


from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import integrator_LangevinBarostat

class LangevinBarostatLocal(integrator_LangevinBarostat):
  def __init__(self, system, rng, temperature):
    'The (local) Velocity Verlet Integrator.'
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, integrator_LangevinBarostat, system, rng, temperature)

if pmi.isController :
  class LangevinBarostat(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.integrator.LangevinBarostatLocal',
      pmiproperty = [ 'gammaP', 'pressure', 'mass' ],
      pmicall = [ "setMassByFrequency" ]
    )
