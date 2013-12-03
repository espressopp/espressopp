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
********************************************************************
**LBInit** - abstract base class for LatticeBoltzmann initialization
********************************************************************

This abstract base class provides the interface and some basic
functionality for classes that (re)initialize populations and handle external forces
  
It provides the following methods:

* createDenVel(rho0,u0)
    sets initial density and velocity 
* setForce()
    sets external force to a specific values
* setForce()
    adds a specific value to the existing forces
"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.LatticeBoltzmann import *
from _espresso import integrator_LBInit

class LBInitLocal(integrator_LBInit):
    """Abstract local base class for LBInit."""
    def createDenVel(self,rho0,u0):
        """createDenVel helps to create initial populations with desired density and
        velocity. By default either a uniform conformation is created by function LBInitPopUniform or 
        a conformation with a constant density and sin-wave-like v_z component as a function of x by
        function LBInitPopWave.
        
        Example:
    
        >>> initPop = espresso.integrator.LBInitPopUniform(system,lb)
        >>> initPop.createDenVel(1.0, Real3D(0.,0.,0.0))
        >>> # first number is the density, second number is a vector of velocity
    
        Example:
        
        >>> initPop = espresso.integrator.LBInitPopWave(system,lb)
        >>> initPop.createDenVel(1.0, Real3D(0.,0.,0.0005))
        >>> # the Real3D vector in this case includes amplitudes of the velocities
        
        """
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.createDenVel(self,rho0,u0)        
    def setForce(self,force):
        """setForce sets an external force onto the system. It is either a constant body force 
        (gravity-like) coded by LBInitConstForce or a sin-wave-like f_z force component as a function
        of x provided by LBInitPeriodicForce.
        
        Example:
    
        >>> lbforce = espresso.integrator.LBInitConstForce(system,lb)
        >>> lbforce.setForce(Real3D(0.,0.,0.0005))
        >>> # a vector sets the external body force directly in lb-units
    
        Example:
        
        >>> lbforce = espresso.integrator.LBInitPeriodicForce(system,lb)
        >>> lbforce.setForce(Real3D(0.,0.,0.0005))
        >>> # a vector sets the external body force amplitude
        
        """
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setForce(self,force)
    def addForce(self,force):
        """addForce adds an external force onto the system. All existing forces will be preserved!
        A user might use it for a superposition of forces desired in a specific application. A constant 
        (gravity-like) force coded by LBInitConstForce and a sin-wave-like f_z force component as a function
        of x provided by LBInitPeriodicForce.
        
        Example:
    
        >>> lbforce = espresso.integrator.LBInitConstForce(system,lb)
        >>> lbforce.addForce(Real3D(0.,0.,0.0005))
        >>> # a vector adds the external body force directly in lb-units
    
        Example:
        
        >>> lbforce = espresso.integrator.LBInitPeriodicForce(system,lb)
        >>> lbforce.addForce(Real3D(0.,0.,0.0005))
        >>> # a vector adds the external body force with a Real3D amplitude
        
        """
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.addForce(self,force)
        
if pmi.isController :
    class LBInit():
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            )
        
