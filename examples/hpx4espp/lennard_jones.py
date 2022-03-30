#!/usr/bin/env python3
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
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--no-hpx4espp", action="store_false", dest="hpx4espp", default=True, help="Use the standard version of modules")
parser.add_argument("--no-hpxStart", action="store_false", dest="hpxStart", default=True, help="Use HPX modules but without multithreading")
args, _ = parser.parse_known_args()

isteps      = 1000
rc          = 2.5
skin        = 0.3
timestep    = 0.005
rho         = 0.8442

# set temperature to None for NVE-simulations
temperature = 1.0

# number of subdomains (should be more than the number of cores/mpi rank)
numSubs     = 32

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
print(espressopp.Version().info())
print('Setting up simulation ...')
particles_per_direction = 32
num_particles = particles_per_direction**3
x, y, z, Lx, Ly, Lz = espressopp.tools.lattice.createCubic(num_particles, rho=rho, perfect=True)
_b, _a, x, y, z, Lx, Ly, Lz = espressopp.tools.replicate([], [], x, y, z, Lx, Ly, Lz, xdim=2, ydim=2, zdim=2)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
box = (Lx, Ly, Lz)

if args.hpx4espp:
    import espressopp.hpx4espp
    system, integrator = espressopp.hpx4espp.standard_system.Default(numSubs=numSubs, box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)
else:
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
if args.hpx4espp:
    vl      = espressopp.hpx4espp.VerletList(system, cutoff = rc)
    potLJ   = espressopp.hpx4espp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
    interLJ = espressopp.hpx4espp.interaction.VerletListLennardJones(vl)
    interLJ.setPotential(type1=0, type2=0, potential=potLJ)
    system.addInteraction(interLJ)
else:
    vl      = espressopp.VerletList(system, cutoff = rc)
    potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
    interLJ = espressopp.interaction.VerletListLennardJones(vl)
    interLJ.setPotential(type1=0, type2=0, potential=potLJ)
    system.addInteraction(interLJ)

# print simulation parameters
print('')
print('number of particles = ', num_particles)
print('density             = ', density)
print('rc                  = ', rc)
print('dt                  = ', integrator.dt)
print('skin                = ', system.skin)
print('temperature         = ', temperature)
print('isteps              = ', isteps)
print('NodeGrid            = ', system.storage.getNodeGrid())
print('CellGrid            = ', system.storage.getCellGrid())
print('')
print('hpx4espp            = ', args.hpx4espp)
print('hpxStart            = ', args.hpxStart)
print('')

if args.hpx4espp:
    hpx = espressopp.hpx4espp.HPXRuntime()

espressopp.tools.analyse.info(system, integrator)

if args.hpx4espp and args.hpxStart:
    hpx.start()

start_time = time.process_time()
integrator.run(isteps)
end_time = time.process_time()

if args.hpx4espp and args.hpxStart:
    hpx.stop()

espressopp.tools.analyse.info(system, integrator)
espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)
