from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import integrator_StochasticVelocityRescaling 

class StochasticVelocityRescalingLocal(integrator_StochasticVelocityRescaling):
    'The (local) StochasticVelocityRescaling Thermostat.'
    def __init__(self, system):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_StochasticVelocityRescaling, system)

if pmi.isController :
    class StochasticVelocityRescaling(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.StochasticVelocityRescalingLocal',
            pmiproperty = [ 'temperature', 'coupling' ]
            )
