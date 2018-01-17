#!/usr/bin/env python2
#
#  Copyright (C) 2013-2017(H)
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
# 
# -*- coding: utf-8 -*-
#
###########################################################################
#                                                                         #
#  This Python script may be used to simulate a monatomic LJ fluid in the #
#  NVE or NVT ensemble. The starting configuration may be taken from      #
#  either a LAMMPS data file or by generating coordinates on a lattice.   #
#                                                                         #
###########################################################################

def openmpi_workaround():

  # find the flag to be set to open all shared libraries at once

  try:
    import dl
    globalFlag = dl.RTLD_GLOBAL
  except:
    try:
      import ctypes
      globalFlag = ctypes.RTLD_GLOBAL
    except:
      print 'ATTENTION: could not find flag RTLD_GLOBAL for dlopen'
      # take a good value for Linux (but is platform-dependent)
      globalFlag = 256

  # now set this flag so that dlopen will use it

  import sys
  flags = sys.getdlopenflags()
  sys.setdlopenflags(flags | globalFlag)

openmpi_workaround()

import sys
import time
import _espressopp
import espressopp
import mpi4py.MPI as MPI
import logging
from espressopp import Real3D, Int3D
from espressopp.tools import lammps
from espressopp.tools import decomp
from espressopp.tools import lattice
from espressopp.tools import timers

# logging.getLogger("Storage").setLevel(logging.INFO)

# simulation parameters (nvt = False implies NVE)
steps = 10
rc = 2.5
skin = 0.3
nvt = False
timestep = 0.005


######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
sys.stdout.write('Setting up simulation ...\n')
x, y, z, Lx, Ly, Lz, vx, vy, vz = lammps.read('espressopp_lennard_jones.start')
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
system = espressopp.System()
system.rng = espressopp.esutil.RNG()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size,size,rc,skin)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# add particles to the system and then decompose
props = ['id', 'type', 'mass', 'pos', 'v']
new_particles = []
for i in range(num_particles):
  part = [i + 1, 0, 1.0, Real3D(x[i], y[i], z[i]), Real3D(vx[i], vy[i], vz[i])]
  new_particles.append(part)
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

# all particles interact via a LJ interaction (use Verlet lists)
vl = espressopp.VerletList(system, cutoff=rc+system.skin)
#potLJ = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=False)
#interLJ = espressopp.interaction.VerletListLennardJones(vl)
potLJ = espressopp.interaction.LennardJonesGromacs(epsilon=1.0, sigma=1.0, r1=2.0, cutoff=rc, shift=False)
potLJX = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=False)
interLJ = espressopp.interaction.VerletListLennardJonesGromacs(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

for i in range(2,101):
  r = (i/100.0) * rc
  print r, potLJ.computeEnergy(r), potLJX.computeEnergy(r),  potLJX.computeForce(Real3D(r, 0, 0))

# setup integrator
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = timestep

if(nvt):
  langevin = espressopp.integrator.LangevinThermostat(system)
  langevin.gamma = 1.0
  langevin.temperature = 1.0
  integrator.addExtension(langevin)

print ''
print 'number of particles =', num_particles
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'nvt =', nvt
print 'steps =', steps
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
temperature = espressopp.analysis.Temperature(system)
pressure = espressopp.analysis.Pressure(system)
pressureTensor = espressopp.analysis.PressureTensor(system)

fmt = '%5d %8.4f %10.5f %8.5f %12.3f %12.3f %12.3f\n'

T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
sys.stdout.write(' step     T          P        Pxy       etotal     epotential    ekinetic\n')
sys.stdout.write(fmt % (0, T, P, Pij[3], Ek + Ep, Ep, Ek))

start_time = time.clock()
integrator.run(steps)
end_time = time.clock()
T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
sys.stdout.write(fmt % (steps, T, P, Pij[3], Ek + Ep, Ep, Ek))
sys.stdout.write('\n')

timers.show(integrator.getTimers(), precision=3)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))
