"""
******************************************************
**LangevinBarostat** - Langevin-Hoover barostat Object
******************************************************

** NOT FINISHED **

This is the barostat implementation to perform Langevin dynamics in a Hoover style extended system
according to the original paper [Quigley04]_.
If LangevinBarostat is defined (as a property of integrator) the integration equations will be 
modified. The volume of system :math:`V` introduced as a dynamical variable:

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

Example:

    >>> rng = espresso.esutil.RNG()
    >>> langevinP = espresso.integrator.LangevinBarostat(system, rng)
    >>> langevinP.gammaP = 0.05
    >>> langevinP.pressure = 1.0
    >>> langevinP.piston = pow(10.0, 4)
    >>> integrator.langevinBarostat = langevinP

**!IMPORTANT** In order to run *npt* simulation one should separately define thermostat as well.

********************************* end

Definition:

    In order to define the Berendsen barostat
    
    >>> berendsen = espresso.integrator.Berendsen(system)
    
    one should have the System_ defined.

.. _System: espresso.System.html

Properties:

*   *berendsen.tau*

    The property 'tau' defines the time parameter :math:`\\tau`.

*   *berendsen.pressure*
    
    The property 'pressure' defines the external pressure :math:`P_{0}`.
    
Setting the integration property:
    
    >>> integrator.berendsen = berendsen
    
    It will define Berendsen barostat as a property of integrator.
    
One more example:

    >>> berendsen_barostat = espresso.integrator.Berendsen(system)
    >>> berendsen_barostat.tau = 10.0
    >>> berendsen_barostat.pressure = 3.5
    >>> integrator.berendsen = berendsen_barostat


Canceling the barostat:
    
    If one do not need the pressure regulation in system anymore or need to switch
    the ensamble or whatever :)

    >>> # define barostat with parameters
    >>> berendsen = espresso.integrator.Berendsen(system)
    >>> berendsen.tau = 0.8
    >>> berendsen.pressure = 15.0
    >>> integrator.berendsen = berendsen
    >>> ...
    >>> # some runs
    >>> ...
    >>> # erase Berendsen barostat
    >>> integrator.berendsen = None
    >>> # the next runs will not include the system size and particle coordinates scaling

References:

.. [Quigley04] D. Quigley, M.I.J. Probert, *J. Chem. Phys.*, 120, **2004**, p. 11432

"""


from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import integrator_LangevinBarostat

class LangevinBarostatLocal(integrator_LangevinBarostat):
  def __init__(self, system, rng):
    'The (local) Velocity Verlet Integrator.'
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, integrator_LangevinBarostat, system, rng)

if pmi.isController :
  class LangevinBarostat(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict( cls =  'espresso.integrator.LangevinBarostatLocal', pmiproperty = [ 'gammaP', 'pressure', 'mass' ] )
