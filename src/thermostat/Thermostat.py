from espresso import pmi

class ThermostatLocal(object):
    pass

if pmi.IS_CONTROLLER:
    class Thermostat(object):
#         __metaclass__ = pmi.Proxy
#         pmiproxydefs = {
#             'pmiproperty': [' temperature' ]
#             }
        pass
        
