#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


"""
****************************************
**espresso.integrator.LatticeBoltzmann**
****************************************

"""
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
