#!/usr/bin/env python

###########################################################################
#                                                                         #
#  This Python script may be used to simulate a monatomic LJ fluid in the #
#  NVE or NVT ensemble. The starting configuration may be taken from      #
#  either a LAMMPS data file or by generating coordinates on a lattice.   #
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
from espresso.tools.init_cfg import lattice

# integration steps, cutoff, skin and thermostat flag (nvt = False is nve)
steps = 1000
rc = 2.5
skin = 0.3
nvt = False

# initial configuration: (1) LAMMPS, (2) lattice or (3) GROMACS
init_cfg = 1

if(init_cfg == 1):
  # LAMMPS with N = 32000
  # useful for checking for identical results against LAMMPS
  x, y, z, Lx, Ly, Lz = lammps.read('data.lj')
  num_particles = len(x)
elif(init_cfg == 2):
  # cubic lattice with user-defined values of N and rho
  # num_particles should be a perfect cube (e.g. 25**3=15625, 32**3=32768)
  num_particles = 20**3
  rho = 0.8442
  x, y, z, Lx, Ly, Lz = lattice.create(num_particles, rho, perfect=False)
else:
  sys.std.write('init_cfg invalid: ' + str(init_cfg) + '. Exiting ...\n')
  sys.exit(1)



######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################

density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
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
integrator.dt = 0.001

if(nvt):
  langevin = espresso.integrator.Langevin(system)
  langevin.gamma = 1.0
  langevin.temperature = 1.0
  integrator.langevin = langevin
  integrator.dt = 0.01

print ''
print 'number of particles =', num_particles
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

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

start_time = time.clock()
integrator.run(steps)
print 'CPU time =', time.clock() - start_time, 's'
T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
sys.stdout.write(fmt % (steps, T, P, Pij[3], Ek + Ep, Ep, Ek))
sys.exit(1)
# comment out line above for production run

intervals = 20
nsteps = steps / intervals
for i in range(1, intervals + 1):
  integrator.run(nsteps)
  step = nsteps * i
  T = temperature.compute()
  P = pressure.compute()
  Pij = pressureTensor.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interLJ.computeEnergy()
  sys.stdout.write(fmt % (step, T, P, Pij[3], Ek + Ep, Ep, Ek))
