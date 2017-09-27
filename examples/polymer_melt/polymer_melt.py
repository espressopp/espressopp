#!/usr/bin/env python2
#  Copyright (C) 2012-2017(H)
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

###########################################################################
#                                                                         #
#  ESPResSo++ Python script for a Polymer Melt System including           #
#  runtime details
#                                                                         #
###########################################################################

import time
import espressopp

nsteps      = 10
isteps      = 100
rc          = pow(2.0, 1.0/6.0)
skin        = 0.4
timestep    = 0.005

# set temperature to None for NVE-simulations
temperature = 1.0

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
print espressopp.Version().info()
print 'Setting up simulation ...'
bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.lammps.read('polymer_melt.lammps')
bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=1, ydim=1, zdim=1)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
box = (Lx, Ly, Lz)
system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)

# add particles to the system and then decompose
# do this in chunks of 1000 particles to speed it up
props = ['id', 'type', 'mass', 'pos']
new_particles = []
for i in range(num_particles):
  part = [i + 1, 0, 1.0, espressopp.Real3D(x[i], y[i], z[i])]
  new_particles.append(part)
  if i % 1000 == 0:
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
    new_particles = []
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

# Lennard-Jones with Verlet list
vl      = espressopp.VerletList(system, cutoff = rc)
potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
interLJ = espressopp.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# FENE bonds
fpl = espressopp.FixedPairList(system.storage)
fpl.addBonds(bonds)
potFENE = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, fpl, potFENE)
system.addInteraction(interFENE)

# Cosine with FixedTriple list
ftl = espressopp.FixedTripleList(system.storage)
ftl.addTriples(angles)
potCosine = espressopp.interaction.Cosine(K=1.5, theta0=3.1415926)
interCosine = espressopp.interaction.FixedTripleListCosine(system, ftl, potCosine)
system.addInteraction(interCosine)

# print simulation parameters
print ''
print 'number of particles = ', num_particles
print 'density             = ', density
print 'rc                  = ', rc
print 'dt                  = ', integrator.dt
print 'skin                = ', system.skin
print 'temperature         = ', temperature
print 'nsteps              = ', nsteps
print 'isteps              = ', isteps
print 'NodeGrid            = ', system.storage.getNodeGrid()
print 'CellGrid            = ', system.storage.getCellGrid()
print ''

# espressopp.tools.decomp.tuneSkin(system, integrator)

espressopp.tools.analyse.info(system, integrator)
start_time = time.clock()
for k in range(nsteps):
  integrator.run(isteps)
  espressopp.tools.analyse.info(system, integrator)
end_time = time.clock()
espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)

