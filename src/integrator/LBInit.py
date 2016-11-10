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


R"""
This abstract class provides the interface to (re-)initialize populations and handle external forces.

.. py:class:: espressopp.integrator.LBInit

    .. py:method:: createDenVel(rho0,u0)
  
        to set initial density and velocity of the LB-fluid.
  
        :param real rho0: density
        :param Real3D u0: velocity
  
    The following options for LB-fluid initialization are supported:
  
    * :class:`espressopp.integrator.LBInitPopUniform` A typical choice. It initializes uniformly distributed density and velocity: On every lattice site the density is ``rho0`` and velocity is ``u0``
    
    * :class:`espressopp.integrator.LBInitPopWave` for uniform density at every lattice site, but harmonic velocity :math:`v_z (x)` with the period of lattice sites in *x*-direction
    
    .. py:method:: setForce(value)

        to set an external force onto LB-fluid.

        :param Real3D value: value of the force

    .. py:method:: addForce(force)

        to add a new external force to the existing one.

        :param Real3D force: value of the force

    Two main external force types are implemented:

    * :class:`espressopp.integrator.LBInitConstForce` to manage constant (gravity-like) forces acting on every lattice site and
    
    * :class:`espressopp.integrator.LBInitPeriodicForce` to manage periodic (sin-like) forces

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
        
