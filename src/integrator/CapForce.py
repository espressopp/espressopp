from espresso.esutil import cxxinit
from espresso import pmi
from espresso.integrator.Extension import *
from _espresso import integrator_CapForce 

class CapForceLocal(ExtensionLocal, integrator_CapForce):
    'The (local) force capping part.'
    def __init__(self, system, capForce, particleGroup = None):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if (particleGroup == None) or (particleGroup.size() == 0):
              cxxinit(self, integrator_CapForce, system, capForce)
            else:
              cxxinit(self, integrator_CapForce, system, capForce, particleGroup)

if pmi.isController :
    class CapForce(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.CapForceLocal',
            pmicall = ['setCapForce', 'setAbsCapForce', 'getCapForce', 'getAbsCapForce'],
            pmiproperty = [ 'particleGroup' ]
            )
