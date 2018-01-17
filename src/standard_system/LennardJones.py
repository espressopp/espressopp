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
***************************************
espressopp.standard_system.LennardJones
***************************************


.. function:: espressopp.standard_system.LennardJones(num_particles, box, rc, skin, dt, epsilon, sigma, shift, temperature, xyzfilename, xyzrfilename)

		:param num_particles: 
		:param box: (default: (000))
		:param rc: (default: 1.12246)
		:param skin: (default: 0.3)
		:param dt: (default: 0.005)
		:param epsilon: (default: 1.0)
		:param sigma: (default: 1.0)
		:param shift: (default: 'auto')
		:param temperature: (default: None)
		:param xyzfilename: (default: None)
		:param xyzrfilename: (default: None)
		:type num_particles: 
		:type box: 
		:type rc: real
		:type skin: real
		:type dt: real
		:type epsilon: real
		:type sigma: real
		:type shift: 
		:type temperature: 
		:type xyzfilename: 
		:type xyzrfilename: 
		
		return random Lennard Jones system and integrator:
		if tempearture is != None then Langevin thermostat is set to temperature (gamma is 1.0)
"""
import espressopp
import mpi4py.MPI as MPI

def LennardJones(num_particles, box=(0,0,0), rc=1.12246, skin=0.3, dt=0.005, epsilon=1.0, sigma=1.0, shift='auto', temperature=None, xyzfilename=None, xyzrfilename=None):

    
  if xyzfilename and xyzrfilename:
     print "ERROR: only one of xyzfilename (only xyz data) or xyzrfilename (additional particle radius data) can be provided."
     sys.exit(1)

  if xyzrfilename: 
    pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf, radiusf = espressopp.tools.readxyzr(xyzrfilename)
    box = (Lxf, Lyf, Lzf)
    num_particles = len(pidf)
  elif xyzfilename:
    pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf = espressopp.tools.readxyz(xyzfilename)
    box = (Lxf, Lyf, Lzf)
    num_particles = len(pidf)
  else:
    if box[0]<=0 or box[1]<=0 or box[2]<=0:
      print "WARNING: no valid box size specified, box size set to (100,100,100) !"    
    
  system         = espressopp.System()
  system.rng     = espressopp.esutil.RNG()
  system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
  system.skin    = skin
  nodeGrid       = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size,box,rc,skin)
  cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
  system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
  interaction    = espressopp.interaction.VerletListLennardJones(espressopp.VerletList(system, cutoff=rc))
  interaction.setPotential(type1=0, type2=0, potential=espressopp.interaction.LennardJones(epsilon, sigma, rc, shift))
  system.addInteraction(interaction)

  integrator     = espressopp.integrator.VelocityVerlet(system)  
  integrator.dt  = dt
  if (temperature != None):
    thermostat             = espressopp.integrator.LangevinThermostat(system)
    thermostat.gamma       = 1.0
    thermostat.temperature = temperature
    integrator.addExtension(thermostat)
  mass = 1.0
  if xyzrfilename: 
    new_particles = []
    props     = ['id', 'type', 'mass', 'pos', 'v', 'radius']
    for idx in xrange(num_particles):
      part = [ pidf[idx], typef[idx], mass,
               espressopp.Real3D(xposf[idx],yposf[idx],zposf[idx]),
               espressopp.Real3D(xvelf[idx],yvelf[idx],zvelf[idx]),
               radiusf[idx] ]
      new_particles.append(part)
      if idx % 1000 == 0:
        system.storage.addParticles(new_particles, *props)
        system.storage.decompose()
        new_particles = []      
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
  elif xyzfilename: 
    new_particles = []
    props     = ['id', 'type', 'mass', 'pos', 'v']
    for idx in xrange(num_particles):
      part = [ pidf[idx], typef[idx], mass,
               espressopp.Real3D(xposf[idx],yposf[idx],zposf[idx]),
               espressopp.Real3D(xvelf[idx],yvelf[idx],zvelf[idx])]
      new_particles.append(part)
      if idx % 1000 == 0:
        system.storage.addParticles(new_particles, *props)
        system.storage.decompose()
        new_particles = []      
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
  else:  
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
