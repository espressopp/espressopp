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
	**********************************************************
	**LiquidGasLB** - class for lattice Boltzmann methods
	**********************************************************
	
	The LiquidGasLB class is an extension to the integrator class of ESPResSo++.
	It creates a simulation box with specified dimensions and allocates necessary
	memory for a lattice Boltzmann simulation. By default we use D3Q19 lattice model
	(in three dimensions and with 19-velocities on the node model).
	
	LiquidGasLB constructor expects 5 parameters (and a system pointer).
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
	
	>>> lb = espressopp.integrator.LiquidGasLB(system, Ni=Int3D(20, 20, 20))
	>>> # creates a cubic box of 20^3 nodes with default spacing parameters in D3Q19 model.
	
	Example
	
	>>> lb = espressopp.integrator.LiquidGasLB(system, Ni=Int3D(30, 20, 20), a = 0.5, tau = 0.5)
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
	
	>>> lb = espressopp.integrator.LiquidGasLB(system, Ni=Int3D(20, 20, 20))
	>>> lb.lbTemp = 0.0000005
	>>> # creates a box of 20^3 nodes with lattice spacing of 1. and timestep of 1. D3Q19 model.
	>>> # then the fluctuations with the temperature of 0.0000005 are initialized.
	
	Example
	
	>>> lb = espressopp.integrator.LiquidGasLB(system, Ni=Int3D(20, 20, 20))
	>>> lb.gamma_b = 0.5
	>>> lb.gamma_s = 0.5
	>>> # creates a box of 20^3 nodes with lattice spacing of 1. and timestep of 1. D3Q19 model.
	>>> # then the bulk and shear gammas are set to 0.5
	
	"""
from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp import Real3D
from espressopp import Int3D
from espressopp.integrator.Extension import *
from _espressopp import integrator_LiquidGasLB

class LiquidGasLBLocal(ExtensionLocal, integrator_LiquidGasLB):
	def __init__(self, system, Ni , a = 1., tau = 1., numDims = 3, numVels = 19):
		if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
			cxxinit(self, integrator_LiquidGasLB, system, Ni, a, tau, numDims, numVels)

if pmi.isController :
	class LiquidGasLB(Extension):
		__metaclass__ = pmi.Proxy
		pmiproxydefs = dict(
												cls =  'espressopp.integrator.LiquidGasLBLocal',
												pmiproperty = [ 'Ni', 'a', 'tau', 'numDims', 'numVels',
																			 'gamma_b', 'gamma_s', 'gamma_odd', 'gamma_even', 'lbTemp']
												)
