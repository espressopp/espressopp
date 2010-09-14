#!/usr/bin/env python

###########################################################################
#                                                                         #
# Example script for LJ simulation with multiple systems                  #
#                                                                         #
###########################################################################

import sys
import time
import espresso
import MPI
import logging
import lammps_file
from espresso import Real3D, Int3D

# read coordinates and box size
x, y, z, Lx, Ly, Lz = lammps_file.read('data.lj')

num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
rc = 2.5
skin = 0.3
comm = MPI.COMM_WORLD
print ""
print "Create two Lennard-Jones systems each with:"
print "  number of particles =", num_particles
print "  density = %.4f" % density
print "  cutoff =", rc

# compute the number of cells on each node
def calcNumberCells(size, nodes, rc):
  ncells = 1
  while size / (ncells * nodes) >= (rc + skin):
    ncells = ncells + 1
  return ncells - 1

nodeGrid = Int3D(1, 1, comm.size)
cellGrid = Int3D(
  calcNumberCells(size[0], nodeGrid[0], rc),
  calcNumberCells(size[1], nodeGrid[1], rc),
  calcNumberCells(size[2], nodeGrid[2], rc)
  )

print ""
print "Both systems are decomposed the same way:"
print "  NodeGrid = %s" % (nodeGrid,)
print "  CellGrid = %s" % (cellGrid,)

# DEFINITION OF SYSTEM1
# see exec() command to save typing when the number of systems is large
system1 = espresso.System()
system1.rng = espresso.esutil.RNG()
system1.bc = espresso.bc.OrthorhombicBC(system1.rng, size)
system1.skin = skin
system1.storage = espresso.storage.DomainDecomposition(system1, comm, nodeGrid, cellGrid)
# add particles to system1
for pid in range(num_particles):
  system1.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
system1.storage.decompose()
# all particles interact via a LJ interaction (use Verlet lists)
vl1 = espresso.VerletList(system1, cutoff=rc+system1.skin)
potLJ1 = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=False)
interLJ1 = espresso.interaction.VerletListLennardJones(vl1)
interLJ1.setPotential(type1=0, type2=0, potential=potLJ1)
system1.addInteraction(interLJ1)
# setup integrator
integrator1 = espresso.integrator.VelocityVerlet(system1)
integrator1.dt = 0.001
# thermostat
langevin1 = espresso.integrator.Langevin(system1)
integrator1.langevin = langevin1
langevin1.gamma = 10.0
langevin1.temperature = 1.0
# analysis
temperature1 = espresso.analysis.Temperature(system1)
pressure1 = espresso.analysis.Pressure(system1)
pressureTensor1 = espresso.analysis.PressureTensor(system1)

# DEFINITION OF SYSTEM2
system2 = espresso.System()
system2.rng = espresso.esutil.RNG()
system2.bc = espresso.bc.OrthorhombicBC(system2.rng, size)
system2.skin = skin
system2.storage = espresso.storage.DomainDecomposition(system2, comm, nodeGrid, cellGrid)
# add particles to system2
for pid in range(num_particles):
  system2.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
system2.storage.decompose()
# all particles interact via a LJ interaction (use Verlet lists)
vl2 = espresso.VerletList(system2, cutoff=rc+system2.skin)
potLJ2 = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=False)
interLJ2 = espresso.interaction.VerletListLennardJones(vl2)
interLJ2.setPotential(type1=0, type2=0, potential=potLJ2)
system2.addInteraction(interLJ2)
# setup integrator
integrator2 = espresso.integrator.VelocityVerlet(system2)
integrator2.dt = 0.001
# thermostat
langevin2 = espresso.integrator.Langevin(system2)
integrator2.langevin = langevin2
langevin2.gamma = 10.0
langevin2.temperature = 2.0
# analysis
temperature2 = espresso.analysis.Temperature(system2)
pressure2 = espresso.analysis.Pressure(system2)
pressureTensor2 = espresso.analysis.PressureTensor(system2)

fmt = '%5d %8.4f %10.5f %8.5f %12.3f %12.3f %12.3f\n'

sys.stdout.write('\n')
sys.stdout.write('Initial temperatures:\n')
sys.stdout.write('  Target T of System 1 = %.1f\n' % langevin1.temperature)
sys.stdout.write('  Target T of System 2 = %.1f\n' % langevin2.temperature)
sys.stdout.write('\n')

T1 = temperature1.compute()
P = pressure1.compute()
Pij = pressureTensor1.compute()
Ek = 0.5 * T1 * (3 * num_particles)
Ep = interLJ1.computeEnergy()
sys.stdout.write(' =============================== SYSTEM 1 ===============================\n')
sys.stdout.write(' step     T        P        Pxy       etotal      epotential     ekinetic\n')
sys.stdout.write(fmt % (0, T1, P, Pij[3], Ek + Ep, Ep, Ek))

numsteps = 20
intervals = 10
for i in range(1, intervals + 1):
  step = i * numsteps
  integrator1.run(numsteps)
  T1 = temperature1.compute()
  P = pressure1.compute()
  Pij = pressureTensor1.compute()
  Ek = 0.5 * T1 * (3 * num_particles)
  Ep = interLJ1.computeEnergy()
  sys.stdout.write(fmt % (step, T1, P, Pij[3], Ek + Ep, Ep, Ek))

T2 = temperature2.compute()
P = pressure2.compute()
Pij = pressureTensor2.compute()
Ek = 0.5 * T2 * (3 * num_particles)
Ep = interLJ2.computeEnergy()
sys.stdout.write('\n')
sys.stdout.write(' =============================== SYSTEM 2 ===============================\n')
sys.stdout.write(' step     T        P        Pxy       etotal      epotential     ekinetic\n')
sys.stdout.write(fmt % (0, T2, P, Pij[3], Ek + Ep, Ep, Ek))

for i in range(1, intervals + 1):
  step = i * numsteps
  integrator2.run(numsteps)
  T2 = temperature2.compute()
  P = pressure2.compute()
  Pij = pressureTensor2.compute()
  Ek = 0.5 * T2 * (3 * num_particles)
  Ep = interLJ2.computeEnergy()
  sys.stdout.write(fmt % (step, T2, P, Pij[3], Ek + Ep, Ep, Ek))

# exchange temperatures between systems
langevin1.temperature = T2
langevin2.temperature = T1

sys.stdout.write('\n')
sys.stdout.write('Exchanges final temperatures:\n')
sys.stdout.write('  Target T of System 1 = %.4f\n' % T2)
sys.stdout.write('  Target T of System 2 = %.4f\n' % T1)

T1 = temperature1.compute()
P = pressure1.compute()
Pij = pressureTensor1.compute()
Ek = 0.5 * T1 * (3 * num_particles)
Ep = interLJ1.computeEnergy()
sys.stdout.write('\n')
sys.stdout.write(' =============================== SYSTEM 1 ===============================\n')
sys.stdout.write(' step     T        P        Pxy       etotal      epotential     ekinetic\n')
sys.stdout.write(fmt % (numsteps*intervals, T1, P, Pij[3], Ek + Ep, Ep, Ek))

for i in range(1, intervals + 1):
  step = i * numsteps + numsteps * intervals
  integrator1.run(numsteps)
  T1 = temperature1.compute()
  P = pressure1.compute()
  Pij = pressureTensor1.compute()
  Ek = 0.5 * T1 * (3 * num_particles)
  Ep = interLJ1.computeEnergy()
  sys.stdout.write(fmt % (step, T1, P, Pij[3], Ek + Ep, Ep, Ek))

T2 = temperature2.compute()
P = pressure2.compute()
Pij = pressureTensor2.compute()
Ek = 0.5 * T2 * (3 * num_particles)
Ep = interLJ2.computeEnergy()
sys.stdout.write('\n')
sys.stdout.write(' =============================== SYSTEM 2 ===============================\n')
sys.stdout.write(' step     T        P        Pxy       etotal      epotential     ekinetic\n')
sys.stdout.write(fmt % (numsteps*intervals, T2, P, Pij[3], Ek + Ep, Ep, Ek))

for i in range(1, 11):
  step = i * numsteps + numsteps * intervals
  integrator2.run(numsteps)
  T2 = temperature2.compute()
  P = pressure2.compute()
  Pij = pressureTensor2.compute()
  Ek = 0.5 * T2 * (3 * num_particles)
  Ep = interLJ2.computeEnergy()
  sys.stdout.write(fmt % (step, T2, P, Pij[3], Ek + Ep, Ep, Ek))
