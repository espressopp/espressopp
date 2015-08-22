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
***********************************************************************************************
**LBInit** - abstract class for LatticeBoltzmann initialization and external force management
***********************************************************************************************

This abstract class provides the interface to (re-)initialize populations and handle external forces.

.. function:: createDenVel(rho0,u0)
  
  to set initial density and velocity of the LB-fluid.
  
  :param rho0: density
  :param u0: velocity
  
  At the moment we support the following options for LB-fluid initialization:
  
  - :class:`espressopp.integrator.LBInitPopUniform` for uniformly distributed density and velocity, i.e. on every lattice site the density is rho0 and velocity is u0;
  - :class:`espressopp.integrator.LBInitPopWave` for uniform density at every lattice site, but harmonic velocity :math:`v_z (x)` with the period of lattice sites in *x*-direction;
    
.. function:: setForce(value)

  to set an external force onto LB-fluid.

  :param value: value of the force
  :type value: Real3D

.. function:: addForce(value)

  to add a new external force to the existing one.

  :param value: value of the force
  :type value: Real3D

Two main external force types are implemented:

- :class:`espressopp.integrator.LBInitConstForce` to manage constant (gravity-like) force acting on every lattice site and
- :class:`espressopp.integrator.LBInitPeriodicForce` to manage harmonic (position-dependent) force



.. function:: espressopp.integrator.LBInit.addForce(force)

		:param force: 
		:type force: 
		:rtype: 

.. function:: espressopp.integrator.LBInit.createDenVel(rho0, u0)

		:param rho0: 
		:param u0: 
		:type rho0: 
		:type u0: 
		:rtype: 

.. function:: espressopp.integrator.LBInit.setForce(force)

		:param force: 
		:type force: 
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.LatticeBoltzmann import *
from _espressopp import integrator_LBInit

class LBInitLocal(integrator_LBInit):
    def createDenVel(self,rho0,u0):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.createDenVel(self,rho0,u0)        
    def setForce(self,force):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setForce(self,force)
    def addForce(self,force):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.addForce(self,force)
        
if pmi.isController :
    class LBInit():
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            )
        
