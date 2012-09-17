from espresso.esutil import cxxinit
from espresso import pmi
from espresso.integrator.Extension import *
from _espresso import integrator_ExtForce 

class ExtForceLocal(ExtensionLocal, integrator_ExtForce):
    'The (local) external force part.'
    def __init__(self, system, extForce, particleGroup = None):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if (particleGroup == None) or (particleGroup.size() == 0):
              cxxinit(self, integrator_ExtForce, system, extForce)
            else:
              cxxinit(self, integrator_ExtForce, system, extForce, particleGroup)

if pmi.isController :
    class ExtForce(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.ExtForceLocal',
            pmicall = ['setExtForce', 'getExtForce'],
            pmiproperty = [ 'particleGroup' ]
            )
