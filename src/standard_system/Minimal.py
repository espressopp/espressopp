#  Copyright (C) 2012,2013,2017(H)
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
espressopp.standard_system.Minimal
**********************************


.. function:: espressopp.standard_system.Minimal(num_particles, box, rc, skin, dt, temperature)

		:param num_particles: 
		:param box: 
		:param rc: (default: 1.12246)
		:param skin: (default: 0.3)
		:param dt: (default: 0.005)
		:param temperature: (default: None)
		:type num_particles: 
		:type box: 
		:type rc: real
		:type skin: real
		:type dt: real
		:type temperature: 
		
		Return minimal system and integrator whithout any interactions defined:
		particles have random positions in box
		if tempearture is != None then Langevin thermostat is set to temperature (gamma is 1.0)
"""
import espressopp
import mpi4py.MPI as MPI

def Minimal(num_particles, box, rc=1.12246, skin=0.3, dt=0.005, temperature=None):

  system         = espressopp.System()
  system.rng     = espressopp.esutil.RNG()
  system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
  system.skin    = skin
  nodeGrid       = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size,box,rc,skin)
  cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
  system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

  integrator     = espressopp.integrator.VelocityVerlet(system)  
  integrator.dt  = dt
  if (temperature != None):
    thermostat             = espressopp.integrator.LangevinThermostat(system)
    thermostat.gamma       = 1.0
    thermostat.temperature = temperature
    integrator.addExtension(thermostat)
  
  props = ['id', 'type', 'mass', 'pos', 'v']
  new_particles = []
  pid = 1
  while pid <= num_particles:
    type = 0
    mass = 1.0
    pos  = system.bc.getRandomPos()
    vel  = espressopp.Real3D(0.0, 0.0, 0.0)
    part = [pid, type, mass, pos, vel]
    new_particles.append(part)
    if pid % 1000 == 0:
        system.storage.addParticles(new_particles, *props)
        system.storage.decompose()
        new_particles = []
    pid += 1
  system.storage.addParticles(new_particles, *props)
  system.storage.decompose()
 
  return system, integrator
