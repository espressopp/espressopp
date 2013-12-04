from espresso.esutil import cxxinit
from espresso import pmi

from _espresso import integrator_FreeEnergyCompensation 

class FreeEnergyCompensationLocal(integrator_FreeEnergyCompensation):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system, center=[]):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_FreeEnergyCompensation, system)
            
            # set center of FreeEnergyCompensation force
            if (center != []):
                self.cxxclass.setCenter(self, center[0], center[1], center[2])

    def addForce(self, itype, filename, type):
            """
            Each processor takes the broadcasted interpolation type,
            filename and particle type
            """
            if pmi.workerIsActive():
                self.cxxclass.addForce(self, itype, filename, type)
                
    def computeCompEnergy(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.computeCompEnergy(self)

if pmi.isController :
    class FreeEnergyCompensation(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.FreeEnergyCompensationLocal',
            pmiproperty = [ 'itype', 'filename'],
            pmicall = ['addForce' , 'computeCompEnergy']
            )
