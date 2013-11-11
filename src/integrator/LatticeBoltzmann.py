from espresso.esutil import cxxinit
from espresso import pmi
from espresso import Real3D
from espresso import Int3D
from espresso.integrator.Extension import *
from _espresso import integrator_LatticeBoltzmann 

class LatticeBoltzmannLocal(ExtensionLocal, integrator_LatticeBoltzmann):
    """The (local) Lattice Boltzmann part.
    
    Creates a simulation box with a specified dimensions and allocates the necessary memory for a 
    lattice Boltzmann simulation. Default values of the parameters include: 
    spacing of the latticea = 1.,
    time spacing = 1.,
    number of dimensions = 3,
    number of velocities at the lattice site = 19. 
    
    Example
    
    >>> lb = espresso.integrator.LatticeBoltzmann(system, Ni=Int3D(20, 20, 20))
    >>> # will create a cubic lattice box of 20 sites with default spacing parameters in D3Q19 model. 
    """
    def __init__(self, system, Ni , a = 1., tau = 1., numDims = 3, numVels = 19):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
	      cxxinit(self, integrator_LatticeBoltzmann, system, Ni, a, tau, numDims, numVels)

if pmi.isController :
    class LatticeBoltzmann(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.LatticeBoltzmannLocal',
            pmiproperty = [ 'Ni', 'a', 'tau', 'numDims', 'numVels', 
                           'gamma_b', 'gamma_s', 'gamma_odd', 'gamma_even', 'lbTemp', 'extForce' ]
            )
