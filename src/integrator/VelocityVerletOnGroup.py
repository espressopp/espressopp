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
*********************************************
**espresso.integrator.VelocityVerletOnGroup**
*********************************************

"""
from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.MDIntegrator import *
from _espresso import integrator_VelocityVerletOnGroup

class VelocityVerletOnGroupLocal(MDIntegratorLocal, integrator_VelocityVerletOnGroup):
    'The (local) Velocity Verlet Integrator.'
    def __init__(self, system, group):
        cxxinit(self, integrator_VelocityVerletOnGroup, system, group)

if pmi.isController :
    class VelocityVerletOnGroup(MDIntegrator):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espresso.integrator.VelocityVerletOnGroupLocal',
            pmiproperty = [ 'langevin' ]
            )
