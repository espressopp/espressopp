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
**********************************
espressopp.standard_system.Default
**********************************


.. py:method:: espressopp.standard_system.Default(box, rc = 1.12246, skin = 0.3, dt = 0.005, temperature = None)

		:param box: 
		:param real rc:
		:param real skin:
		:param real dt:
		:param temperature:
		:type box:
		:type temperature: 
		 
		Return default system and integrator, no interactions, no particles are set
		if tempearture is != None then Langevin thermostat is set to temperature (gamma is 1.0)
"""
import espressopp
import mpi4py.MPI as MPI

def Default(box, rc=1.12246, skin=0.3, dt=0.005, temperature=None):

  system         = espressopp.System()
  system.rng     = espressopp.esutil.RNG()
  system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
  system.skin    = skin
  nodeGrid       = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
  cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
  system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

  print "nodeGrid: ",nodeGrid, " cellGrid: ",cellGrid

  integrator     = espressopp.integrator.VelocityVerlet(system)  
  integrator.dt  = dt
  if (temperature != None):
    thermostat             = espressopp.integrator.LangevinThermostat(system)
    thermostat.gamma       = 1.0
    thermostat.temperature = temperature
    integrator.addExtension(thermostat)
   
  return system, integrator
