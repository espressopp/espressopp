#!/usr/bin/env python
# -*- coding: utf-8 -*-



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
from espresso.tools import timers

# integration steps, cutoff, skin and thermostat flag (nvt = False is nve)
steps = 1000
rc = 2.5
skin = 0.3
nvt = False
timestep = 0.005

# initial configuration: (1) LAMMPS, (2) lattice or (3) GROMACS
init_cfg = 2

if(init_cfg == 1):
  # LAMMPS with N = 32000
  # useful for checking for identical results against LAMMPS
  file = sys.path[0][:sys.path[0].find('espressopp')] + 'espressopp/examples/data.lj'
  x, y, z, Lx, Ly, Lz = lammps.read(file)
  num_particles = len(x)
elif(init_cfg == 2):
  # cubic lattice with user-defined values of N and rho
  # num_particles should be a perfect cube (e.g. 25**3=15625, 32**3=32768)
  num_particles = 20**3
  rho = 0.8442 # density
  x, y, z, Lx, Ly, Lz = lattice.create(num_particles, rho, perfect=False)
else:
  sys.stdout.write('init_cfg invalid: ' + str(init_cfg) + '. Exiting ...\n')
  sys.exit(1)



######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
sys.stdout.write('Setting up simulation ...\n')
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

# compute the number of cells on each node
def calcNumberCells(size, nodes, cutoff):
    ncells = 1
    while size / (ncells * nodes) >= cutoff:
       ncells = ncells + 1
    return ncells - 1

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
#nodeGrid = Int3D(1, 1, comm.size)
#cellGrid = Int3D(
    #calcNumberCells(size[0], nodeGrid[0], rc),
    #calcNumberCells(size[1], nodeGrid[1], rc),
    #calcNumberCells(size[2], nodeGrid[2], rc)
    #)


system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)



# add particles to the system and then decompose
for pid in range(num_particles):
  system.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
system.storage.decompose()



# all particles interact via a LJ interaction (use Verlet lists)
vl = espresso.VerletList(system, cutoff=rc+system.skin)

potLJ = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=False)
#potLJ = espresso.interaction.SoftCosine(A=1.0, cutoff=rc, shift=False)

interLJ = espresso.interaction.VerletListLennardJones(vl)
#interLJ = espresso.interaction.VerletListSoftCosine(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)

system.addInteraction(interLJ)


# setup integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = timestep

if(nvt):
  langevin = espresso.integrator.Langevin(system)
  langevin.gamma = 1.0
  langevin.temperature = 1.0
  integrator.langevin = langevin
  integrator.dt = timestep

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
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)
pressureTensor = espresso.analysis.PressureTensor(system)

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
