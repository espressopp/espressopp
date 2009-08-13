from espresso import pmi
from espresso.esutil import *
from espresso.thermostat import *
from _espresso import thermostat_Langevin

class LangevinLocal(ThermostatLocal, thermostat_Langevin):
    def __init__(self, temperature, gamma):
        cxxinit(self, thermostat_Langevin, temperature, gamma)

if pmi.IS_CONTROLLER:
    class Langevin(Thermostat):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = \
            dict(cls='espresso.thermostat.LangevinLocal',
                 pmicall = ['connect', 'disconnect'],
                 pmiproperty = [ 'gamma' ])
