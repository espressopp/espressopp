#!/bin/python
# -*- coding: iso-8859-1 -*-
import espresso
import MPI
from espresso.tools import decomp
from espresso import Real3D
# tes# testt
size=(20,20,20)

rc = 2.5
skin = 0.3
num_particles = 5
timestep = 1

system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espresso.storage.DomainDecomposition(system, comm, nodeGrid, cellGrid)
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)

props = ['id', 'pos', 'v']
parts = []
for pid in range(num_particles):
  parts.append([pid, Real3D(0, 0, 2.5*pid), (0,0,1) ])

system.storage.addParticles(parts, *props)

grp=espresso.ParticleGroup(system.storage)

grp.add(1)
grp.add(2)
grp.add(num_particles - 1)

grp.show()

integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = timestep

langevin = espresso.integrator.Langevin(system)
langevin.gamma = 1.0
langevin.temperature = 1.0
integrator.langevin = langevin
integrator.dt = timestep
  
configurations = espresso.analysis.Configurations(system)
configurations.gather()

for i in range(1,10):
  configurations.clear()
  integrator.run(200)
  configurations.gather()

  for conf in configurations:
     NP = conf.size
     print("Configuration : %d particles\n"%NP)
     for id in conf:
        pos = conf[id]
        print("%d : %g %g %g\n"%(id, pos.x, pos.y, pos.z))

  grp.show()
