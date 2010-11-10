#!/usr/bin/env python

###########################################################################
#                                                                         #
#  Example script for LJ simulation with or without Langevin thermostat   #
#                                                                         #
###########################################################################

import sys
import time
import espresso
import MPI
import logging
from espresso import Real3D, Int3D
from espresso.tools.convert import lammps
from espresso.tools import decomp

# benchmark or production (bench = True is a short job)
bench = True

steps = 1000
x, y, z, Lx, Ly, Lz = lammps.read('data.lj')
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
rc = 2.5
skin = 0.3

print 'number of particles =', num_particles
print 'density =', density

system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espresso.storage.DomainDecomposition(system, comm, nodeGrid, cellGrid)

# add particles to the system
for pid in range(num_particles):
  system.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))

# assign particles to processors then cells
system.storage.decompose()

# all particles interact via a LJ interaction (use Verlet lists)
vl = espresso.VerletList(system, cutoff=rc+system.skin)
potLJ = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=False)
interLJ = espresso.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# setup integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = 0.005

# analysis
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)
pressureTensor = espresso.analysis.PressureTensor(system)

fmt = '%5d %8.4f %10.5f %8.5f %12.3f %12.3f %12.3f\n'

T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
sys.stdout.write(' step     T        P        Pxy       etotal     epotential    ekinetic\n')
sys.stdout.write(fmt % (0, T, P, Pij[3], Ek + Ep, Ep, Ek))

time1 = time.time()
integrator.run(steps)
time2 = time.time()
T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
sys.stdout.write(fmt % (steps, T, P, Pij[3], Ek + Ep, Ep, Ek))
print("walltime = %g"%(time2 - time1))
