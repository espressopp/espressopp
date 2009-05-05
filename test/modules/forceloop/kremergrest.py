#  Python example script that does the same as forceloop but with
#  epxorted routines from the C++ classes
#
#  Authors: Jon + Thomas
#
#  Done:   
#     - bindings for particle/All and particel/Set
#     - addtional bindings for integrator/VerlocityVerlet 
#     - ParticleWriter as PythonComputer 
#
#  Furter ToDo's 
#     - Extend all C++ classes in xxx/__init__.py as Local + PMI classes
#     - writing UnitTest in Python

import espresso

from _espresso import particles_Storage as Storage
from _espresso import Real3DProperty
from _espresso import Real3D
from _espresso import bc_PBC as PBC
from _espresso import particles_PythonComputer 
from _espresso import particles_All
from _espresso import pairs_All as All2
from _espresso import pairs_List as List
from _espresso import interaction_LennardJones as LennardJones
from _espresso import interaction_FENE as FENE
from _espresso import integrator_VelocityVerlet as VelocityVerlet

##  Simple class to write position of particles

class ParticleWriter(particles_PythonComputer):

    def __init__(self, _property, storage):
        particles_PythonComputer.__init__(self, storage)
        self.property = _property

    def each(self, pid):
        print("%d %s" % (pid, self.property[pid]))

import random

particleStorage = Storage()
position = Real3DProperty(particleStorage)
velocity = Real3DProperty(particleStorage)
force    = Real3DProperty(particleStorage)

N = 3
SIZE = 4.0

for i in range(N):
   for j in range(N):
      for k in range(N):
     
       r = 0.4 + 0.2 * random.random()
       x = (i + r) / N * SIZE;
       y = (j + r) / N * SIZE; 
       z = (k + r) / N * SIZE;

       id = particleStorage.addParticle()
       print 'id = ', id

       position[id] = Real3D(x, y, z);
       velocity[id] = Real3D(x, y, z);
       force[id] = Real3D(0.0)

#  Write all the particles of the storage

#  We could print the particles already here
#  particleStorage.foreach(ParticleWriter(position, particleStorage))

pbc = PBC()
pbc.set(SIZE)

allSet = particles_All(particleStorage)
allSet.foreach(ParticleWriter(position, particleStorage))

# Pair List for nonbonded interactions

allPairs = All2(pbc, allSet, position)

# Bond List for bonded interactions

bondList = List(pbc, particleStorage, position)

bondList.addPair(1,2)
bondList.addPair(2,3)

bondList.addPair(4,5)
bondList.addPair(5,6)

ljint = LennardJones()
ljint.set(1.0, 1.0, 2.5)

fene = FENE()
fene.set(1.0, 0.5, 0.1)
# fene = FENE(r0 = 0.5, K=1.0, rMax = 0.1)

# allPairs.foreach(forcecompute)

integrator = VelocityVerlet(allSet, position, velocity, force)

integrator.setTimeStep(0.005)

integrator.addForce(ljint, allPairs)
integrator.addForce(fene, bondList)

integrator.run(100)

allSet.foreach(ParticleWriter(force, particleStorage))

