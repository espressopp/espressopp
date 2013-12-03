"""
*****************************************************
**BerendsenThermostat** - Berendsen thermostat Object
*****************************************************

This is the Berendsen thermostat implementation according to the original paper [Berendsen84]_.
If Berendsen thermostat is defined (as a property of integrator) then at the each run the system size
and the particle coordinates will be scaled by scaling parameter :math:`\lambda` according to
the formula:

.. math::   \lambda = [1 + \Delta t/\\tau_{T} (T_{0}/T - 1)]^{1/2}

where :math:`\Delta t` - integration timestep, :math:`\\tau_{T}` - time parameter (coupling parameter),
:math:`T_{0}` - external temperature and :math:`T` - instantaneous temperature. 

Example:

    >>> berendsenT = espresso.integrator.BerendsenThermostat(system)
    >>> berendsenT.tau = 1.0
    >>> berendsenT.temperature = 1.0
    >>> integrator.addExtension(berendsenT)


Definition:

    In order to define the Berendsen thermostat
    
    >>> berendsenT = espresso.integrator.BerendsenThermostat(system)
    
    one should have the System_ defined.

.. _System: espresso.System.html

Properties:

*   *berendsenT.tau*

    The property 'tau' defines the time parameter :math:`\\tau_{T}`.

*   *berendsenT.temperature*
    
    The property 'temperature' defines the external temperature :math:`T_{0}`.
    
Setting the integration property:
    
    >>> integrator.addExtension(berendsenT)
    
    It will define Berendsen thermostat as a property of integrator.
    
One more example:

    >>> berendsen_thermostat = espresso.integrator.BerendsenThermostat(system)
    >>> berendsen_thermostat.tau = 0.1
    >>> berendsen_thermostat.temperature = 3.2
    >>> integrator.addExtension(berendsen_thermostat)


Canceling the thermostat:
    
    >>> # define thermostat with parameters
    >>> berendsen = espresso.integrator.BerendsenThermostat(system)
    >>> berendsen.tau = 2.0
    >>> berendsen.temperature = 5.0
    >>> integrator.addExtension(berendsen)
    >>> ...
    >>> # some runs
    >>> ...
    >>> # disconnect Berendsen thermostat
    >>> berendsen.disconnect()

    
  Connecting the thermostat back after the disconnection
  
    >>> berendsen.connect()

"""


from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *
from _espresso import integrator_BerendsenThermostat

class BerendsenThermostatLocal(ExtensionLocal, integrator_BerendsenThermostat):
  def __init__(self, system):
    'The (local) Velocity Verlet Integrator.'
    if not (pmi._PMIComm and pmi._PMIComm.isActive()) or \
            pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      cxxinit(self, integrator_BerendsenThermostat, system)

if pmi.isController:
  class BerendsenThermostat(Extension):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espresso.integrator.BerendsenThermostatLocal',
      pmiproperty = [ 'tau', 'temperature' ]
    )
