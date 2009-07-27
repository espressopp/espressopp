from espresso import pmi
from espresso.esutil import *
from espresso.thermostat import *
from _espresso import thermostat_Langevin

class LangevinLocal(ThermostatLocal, thermostat_Langevin):
    def __init__(self, temperature, gamma):
        cxxinit(self, thermostat_Langevin, temperature, gamma)

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.thermostat')
    class Langevin(Thermostat):
#         __metaclass__ = pmi.Proxy
#         pmiproxydefs = Thermostat.pmiproxydefs
# #        pmiproxydefs['class'] = 'espresso.thermostat.LangevinLocal'
#         pmiproxydefs['pmicall'] = ['connect', 'disconnect']
#         pmiproxydefs['pmiproperty'].extend(['gamma'])

        def __init__(self, temperature, gamma):
            pmiinit(self, 'espresso.thermostat.LangevinLocal',
                    temperature, gamma)

        def connect(self, integrator):
            pmi.call(self.pmiobject.connect, integrator.pmiobject)
