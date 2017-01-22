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



r"""
******************************************
espressopp.interaction.LennardJonesGromacs
******************************************
         
if :math:`d^2 > r_1^2`

.. math::        

	U = 4 \varepsilon (\frac{sigma^{12}}{d^{12}} - \frac{sigma^6}{d^6}) + (d-r_1)^3 (ljsw3 + ljsw4 (d-r_1) + ljsw5)      	

else

.. math::

	U = 4 \varepsilon (\frac{\sigma^{12}}{d^{12}} - \frac{\sigma^6}{d^6})






.. function:: espressopp.interaction.LennardJonesGromacs(epsilon, sigma, r1, cutoff, shift)

		:param epsilon: (default: 1.0)
		:param sigma: (default: 1.0)
		:param r1: (default: 0.0)
		:param cutoff: (default: infinity)
		:param shift: (default: "auto")
		:type epsilon: real
		:type sigma: real
		:type r1: real
		:type cutoff: 
		:type shift: 

.. function:: espressopp.interaction.VerletListLennardJonesGromacs(vl)

		:param vl: 
		:type vl: 

.. function:: espressopp.interaction.VerletListLennardJonesGromacs.getPotential(type1, type2)

		:param type1: 
		:param type2: 
		:type type1: 
		:type type2: 
		:rtype: 

.. function:: espressopp.interaction.VerletListLennardJonesGromacs.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.CellListLennardJonesGromacs(stor)

		:param stor: 
		:type stor: 

.. function:: espressopp.interaction.CellListLennardJonesGromacs.setPotential(type1, type2, potential)

		:param type1: 
		:param type2: 
		:param potential: 
		:type type1: 
		:type type2: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLennardJonesGromacs(system, vl, potential)

		:param system: 
		:param vl: 
		:param potential: 
		:type system: 
		:type vl: 
		:type potential: 

.. function:: espressopp.interaction.FixedPairListLennardJonesGromacs.setPotential(potential)

		:param potential: 
		:type potential: 
"""


"""
    real sig2 = sigma * sigma;
    real sig6 = sig2 * sig2 * sig2;
    ff1 = 48 \varepsilon \sigma^{12}
    ff2 = 24 \varepsilon \sigma^6
    ef1 =  4 \varepsilon \sigma^{12}
    ef2 =  4 \varepsilon \sigma^6
    r1sq = r_^2
    real t = r_c - r_1
	real r6inv = \frac{1}{r_c^6}
	real r8inv = \frac{1}{r_c^8}
	real t2inv = \frac{1}{(r_c - r_1)^2}
	real t3inv = \frac{1}{(r_c - r_1)^3}
	real t3 = (r_c - r_1)^3
	real a6 = \frac{7 r_1 - 10 r_c}{(r_c - r_1)^2 r_c^8}
	real b6 = \frac{9 r_c - 7 r_1}{(r_c - r_1)^3 r_c^8};
	real a12 = \frac{13 r_1 - 16 r_c}{(r_c - r_1)^2 r_c^{14}}
	real b12 = \frac{15 r_c - 13 r_1}{(r_c - r_1)^3 r_c^{14}}
	real c6 = \frac{1}{r_c^6} - (r_c - r_1)^3(\frac{42 r_1 - 60 r_c}{3(r_c - r_1)^2 r_c^8} + \frac{(54 r_c - 42 r_1)(r_c - r_1)}{4(r_c - r_1)^3 r_c^8});
	real c12 = \frac{1}{r_c^{12}} - (r_c - r_1)^3(\frac{156 r_1 - 192 r_c}{3(r_c - r_1)^2 r_c^{14}} + \frac{(180 r_c - 156 r_1)(r_c - r_1}{4(r_c - r_1)^3 r_c^{14}});



	ljsw3 = -4 \varepsilon \sigma^{12} (\frac{156 r_1 - 192 r_c}{3(r_c - r_1)^2 r_c^{14}}) + 4 \varepsilon \sigma^6 \frac{42 r_1 - 60 r_c}{3(r_c - r_1)^2 r_c^8}
	ljsw4 = -4 \varepsilon \sigma^{12} (\frac{180 r_c - 156 r_1}{4(r_c - r_1)^3 r_c^{14}}) + 4 \varepsilon \sigma^6 \frac{54 r_c - 42 r_1}{4(r_c - r_1)^3 r_c^8}
	
	
	
	ljsw5 = -4 \varepsilon \sigma^{12} (\frac{1}{r_c^{12}} - (r_c - r_1)^3(\frac{156 r_1 - 192 r_c}{3(r_c - r_1)^2 r_c^{14}} + \frac{(180 r_c - 156 r_1)(r_c - r_1}{4(r_c - r_1)^3 r_c^{14}})) + 4 \varepsilon \sigma^6 \frac{1}{r_c^6} - (r_c - r_1)^3(\frac{42 r_1 - 60 r_c}{3(r_c - r_1)^2 r_c^8} + \frac{(54 r_c - 42 r_1)(r_c - r_1)}{4(r_c - r_1)^3 r_c^8})


	U = 4 \varepsilon (\frac{\sigma^{12}}{d^{12}} - \frac{\sigma^6}{d^6}) + (d-r_1)^3 ((((-4 \varepsilon \sigma^{12} (\frac{156 r_1 - 192 r_c}{3(r_c - r_1)^2 r_c^{14}}) + 4 \varepsilon \sigma^6 \frac{42 r_1 - 60 r_c}{3(r_c - r_1)^2 r_c^8}
	) + (-4 \varepsilon \sigma^{12} (\frac{180 r_c - 156 r_1}{4(r_c - r_1)^3 r_c^{14}}) + 4 \varepsilon \sigma^6 \frac{54 r_c - 42 r_1}{4(r_c - r_1)^3 r_c^8}
	) (d-r_1) + (-4 \varepsilon \sigma^{12} (\frac{1}{r_c^{12}} - (r_c - r_1)^3(\frac{156 r_1 - 192 r_c}{3(r_c - r_1)^2 r_c^{14}} + \frac{(180 r_c - 156 r_1)(r_c - r_1}{4(r_c - r_1)^3 r_c^{14}}))) + 4 \varepsilon \sigma^6 \frac{1}{r_c^6} - (r_c - r_1)^3(\frac{42 r_1 - 60 r_c}{3(r_c - r_1)^2 r_c^8} + \frac{(54 r_c - 42 r_1)(r_c - r_1)}{4(r_c - r_1)^3 r_c^8})
	)))       
"""
from espressopp import pmi, infinity
from espressopp.esutil import *

from espressopp.interaction.Potential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_LennardJonesGromacs, \
                      interaction_VerletListLennardJonesGromacs, \
                      interaction_CellListLennardJonesGromacs, \
                      interaction_FixedPairListLennardJonesGromacs

class LennardJonesGromacsLocal(PotentialLocal, interaction_LennardJonesGromacs):

    def __init__(self, epsilon=1.0, sigma=1.0, r1=0.0,
                 cutoff=infinity, shift="auto"):
        """Initialize the local LennardJonesGromacs object."""
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if shift =="auto":
                cxxinit(self, interaction_LennardJonesGromacs,
                        epsilon, sigma, r1, cutoff)
            else:
                cxxinit(self, interaction_LennardJonesGromacs,
                        epsilon, sigma, r1, cutoff, shift)

class VerletListLennardJonesGromacsLocal(InteractionLocal, interaction_VerletListLennardJonesGromacs):

    def __init__(self, vl):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_VerletListLennardJonesGromacs, vl)

    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

    def getPotential(self, type1, type2):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getPotential(self, type1, type2)

class CellListLennardJonesGromacsLocal(InteractionLocal, interaction_CellListLennardJonesGromacs):

    def __init__(self, stor):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_CellListLennardJonesGromacs, stor)
        
    def setPotential(self, type1, type2, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, type1, type2, potential)

class FixedPairListLennardJonesGromacsLocal(InteractionLocal, interaction_FixedPairListLennardJonesGromacs):

    def __init__(self, system, vl, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_FixedPairListLennardJonesGromacs, system, vl, potential)
        
    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class LennardJonesGromacs(Potential):
        'The LennardJonesGromacs potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.LennardJonesGromacsLocal',
            pmiproperty = ['epsilon', 'sigma', 'r1']
            )

    class VerletListLennardJonesGromacs(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.VerletListLennardJonesGromacsLocal',
            pmicall = ['setPotential','getPotential']
            )
    class CellListLennardJonesGromacs(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.CellListLennardJonesGromacsLocal',
            pmicall = ['setPotential']
            )
    class FixedPairListLennardJonesGromacs(Interaction):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.interaction.FixedPairListLennardJonesGromacsLocal',
            pmicall = ['setPotential']
            )
