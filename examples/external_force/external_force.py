#!/usr/bin/env python2
#  Copyright (C) 2015-2017(H)
#      Max Planck Institute for Polymer Research
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

#####################################################################################
#                                                                                   #
#  ESPResSo++ Python script for adding an external force to a L-J standard system   #
#                                                                                   #
#####################################################################################

import espressopp
import math

# specify number of particles in cubic box
N = 4
num_particles = pow(N, 3)
# create particle positions in 3d lattice
x, y, z, Lx, Ly, Lz = espressopp.tools.lattice.createCubic(N=num_particles, rho=0.5)

# setup LennardJones system (with 0 particles)
system, integrator = espressopp.standard_system.LennardJones(0, box=(4*Lx, Ly, Lz))

#add particles to the system
particles     = []
property_list = ['id', 'pos']
for k in range(num_particles):
  particles.append([k, espressopp.Real3D(2*Lx+x[k], y[k], z[k])])
  # add particles in chunks of 1000 (this is faster than doing it one by one)
  if k % 1000 == 0:
    system.storage.addParticles(particles, *property_list)
    system.storage.decompose()
    particles = []
system.storage.addParticles(particles, *property_list)
system.storage.decompose()

# specify particles that will be affected by external field
# if external field should be applied to all particles of the system set ef_particleGroup=None
ef_particleGroup  = espressopp.ParticleGroup(system.storage)
for k in range(num_particles):
  particles.append([k, espressopp.Real3D(x[k], y[k], z[k])])
  if k%N>0 and k%N<N-1 :
    ef_particleGroup.add(k)

# define external force
ext_force = espressopp.integrator.ExtForce(system, espressopp.Real3D(1,0,0), ef_particleGroup)

# add external force to the integrator
integrator.addExtension(ext_force)

sock = espressopp.tools.vmd.connect(system)
for i in range(1000):
  # make 10 Velocity-Verlet integration steps
  integrator.run(10)
  # update postions in VMD
  espressopp.tools.vmd.imd_positions(system, sock)
  # modify external force field
  ef = 5.0 * math.cos(0.01*i)
  ext_force.setExtForce(espressopp.Real3D(ef, 0,0))

