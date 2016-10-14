#  Copyright (C) 2012-2016
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008-2011
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


r"""
************************************************************
**LatticeBoltzmann** - class for lattice Boltzmann methods
************************************************************

The LatticeBoltzmann (LB) class is an extension to the integrator class of
ESPResSo++. The main purpose of the LB-fluid in our simulation package is NOT in
fluid dynamics applications or investigation of fluid-solid interfacial
phenomena. We aim at complex soft matter systems, where the LB-fluid is a bulk
solvent and therefore one has rather use some MD particles as solutes. Examples
of such systems range from colloids (point-like MD-particles) to polymer chains
(point-like MD-particles connected into chains) dissolved in some solvent
(LB-fluid) with specific static and dynamic properties.

It is therefore done ON PURPOSE that the user specifies parameters for LB-fluid
in Lennard-Jones (LJ) units. In the kernel of the C++ code we transform these
into LB-units, if neccessary. Such strategy helps users coming from
MD-background to think of the LB-fluid as if it has particle-based structure: to
mimic the solvent one only has to specify such parameters as liquid density,
:math:`rho`, temperature, :math:`T`, and viscosity, :math:`\eta`. For a standard
LJ-fluid one has: :math:`rho \sim 1 [\sigma^{-3}]`, :math:`T \sim 1 [\epsilon]`,
and :math:`\eta \sim 5 [units]`.

.. note::

  Experienced LB-users may find our approach unusual. However, we kindly ask
  them for a feedback, as for us it is also quite novel. Particularly, we are
  interested in suggestions on expansion of the LB-possibilities and would like
  at first get an overview of "what do the people need?". Being it either
  BGK-scheme, implementation of boundary conditions or something else.

It creates a simulation box with specified dimensions and allocates necessary 
memory for a lattice Boltzmann simulation. By default we use D3Q19 lattice model 
(in three dimensions and with 19-velocities on the node model).

LatticeBoltzmann constructor expects 5 parameters (and a system pointer).
These are: lattice size in 3D Ni, lattice spacing a, lattice timestep tau,
number of dimensions and number of velocity vectors on a lattice node.
The lattice size, Ni, is an obligatory parameter and must be set at the
beginning of the simulation.

The default lattice model is D3Q19 (numDims = 3, numVels = 19) and both lattice 
spacing and timestep are set to 1.

Note that at the present stage of development we aim at D3Q19 model.
If you want to use something else, please, feel free to modify the code.

Originally, we had planned this module to operate in 3D only, so if you
need a 2D version, there is a bit more tuning involved. On the other hand,
adding different 3D lattice models (such as D3Q15 or D3Q27) is rather
straightforward.

Example

>>> lb = espressopp.integrator.LatticeBoltzmann(system, Ni=Int3D(20, 20, 20))
>>> # creates a cubic box of 20^3 nodes with default spacing parameters in D3Q19 model.

Example

>>> lb = espressopp.integrator.LatticeBoltzmann(system, Ni=Int3D(30, 20, 20), a = 0.5, tau = 0.5)
>>> # creates a box of 30*20*20 nodes with lattice spacing of 0.5 and timestep of 0.5.
>>> # The model of the lattice is D3Q19.


After initialization of the Lattice Boltzmann module, one has a possibility to
set several properties of the system:

gamma_b and gamma_s are bulk and shear gammas (default values are 0.);

gamma_odd and gamma_even are (hey-hey, surprise!) odd and even gammas (defaults 0.);

(if you are unsure what these gammas are, please refer to any lattice Boltzmann review.
In short, they control correspondent viscosities of the liquid.)

lbTemp is the temperature in lb units for setting up fluctuations (default is 0.);

Example

>>> lb = espressopp.integrator.LatticeBoltzmann(system, Ni=Int3D(20, 20, 20))
>>> lb.lbTemp = 0.0000005
>>> # creates a box of 20^3 nodes with lattice spacing of 1. and timestep of 1. D3Q19 model.
>>> # then the fluctuations with the temperature of 0.0000005 are initialized.

Example

>>> lb = espressopp.integrator.LatticeBoltzmann(system, Ni=Int3D(20, 20, 20))
>>> lb.gamma_b = 0.5
>>> lb.gamma_s = 0.5
>>> # creates a box of 20^3 nodes with lattice spacing of 1. and timestep of 1. D3Q19 model.
>>> # then the bulk and shear gammas are set to 0.5


.. function:: espressopp.integrator.LatticeBoltzmann(system, nodeGrid, Ni, a, tau, numDims, numVels)

		:param system: 
		:param nodeGrid: 
		:param Ni: 
		:param a: (default: 1.)
		:param tau: (default: 1.)
		:param numDims: (default: 3)
		:param numVels: (default: 19)
		:type system: 
		:type nodeGrid: 
		:type Ni: 
		:type a: 
		:type tau: 
		:type numDims: int
		:type numVels: int
"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp import Real3D
from espressopp import Int3D
from espressopp.integrator.Extension import *
from _espressopp import integrator_LatticeBoltzmann 

class LatticeBoltzmannLocal(ExtensionLocal, integrator_LatticeBoltzmann):
    def __init__(self, system, nodeGrid, a = 1., tau = 1., numDims = 3, numVels = 19):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
	      cxxinit(self, integrator_LatticeBoltzmann, system, nodeGrid, a, tau, numDims, numVels)

if pmi.isController :
	class LatticeBoltzmann(Extension):
		__metaclass__ = pmi.Proxy
		pmiproxydefs = dict(
												cls =  'espressopp.integrator.LatticeBoltzmannLocal',
												pmiproperty = ['nodeGrid', 'a', 'tau', 'numDims', 'numVels',
																			 'visc_b','visc_s','gamma_b', 'gamma_s', 'gamma_odd', 'gamma_even',
																			 'lbTemp', 'fricCoeff', 'nSteps', 'profStep', 'getMyNi'],
                                    pmicall = ["getLBMom","setLBMom","saveLBConf","readLBConf","keepLBDump"]
            )
