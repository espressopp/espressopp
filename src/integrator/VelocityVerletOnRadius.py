"""
**********************************************
**espresso.integrator.VelocityVerletOnRadius**
**********************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi
from espresso.integrator.Extension import *
from _espresso import integrator_VelocityVerletOnRadius 

class VelocityVerletOnRadiusLocal(ExtensionLocal, integrator_VelocityVerletOnRadius):
    'The (local) VelocityVerletOnRadius.'
    def __init__(self, system, dampingmass):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_VelocityVerletOnRadius, system, dampingmass)

if pmi.isController :
    class VelocityVerletOnRadius(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.VelocityVerletOnRadiusLocal',
            pmiproperty = [ 'radialDampingMass' ]
            )
