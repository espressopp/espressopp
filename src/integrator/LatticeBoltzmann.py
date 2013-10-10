from espresso.esutil import cxxinit
from espresso import pmi
from espresso import Real3D
from espresso import Int3D
from espresso.integrator.Extension import *
from _espresso import integrator_LatticeBoltzmann 

class LatticeBoltzmannLocal(ExtensionLocal, integrator_LatticeBoltzmann):
    'The (local) Lattice Boltzmann part.'
    def __init__(self, system, Ni , a = 1., tau = 1., rho0 = 1., u0 = Real3D(0.,0.,0.), numDims = 3, numVels = 19):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
#	    if (numDims != None) and (numVels != None):
	      cxxinit(self, integrator_LatticeBoltzmann, system, Ni, a, tau, rho0, u0, numDims, numVels)
#	    elif (numDims == None) and (numVels == None) and (a != None) and (tau != None) and (rho0 != None) and (u0 != None):
#	      cxxinit(self, integrator_LatticeBoltzmann, system, Nx, Ny, Nz, a, tau, rho0, u0)
#	    else:
#	      cxxinit(self, integrator_LatticeBoltzmann, system, Nx, Ny, Nz)    

if pmi.isController :
    class LatticeBoltzmann(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LatticeBoltzmannLocal',
            pmiproperty = [ 'Ni', 'a', 'tau', 'rho0', 'u0', 'numDims', 'numVels', 
                           'gamma_b', 'gamma_s', 'gamma_odd', 'gamma_even', 'lbTemp', 'extForce' ]
            )
