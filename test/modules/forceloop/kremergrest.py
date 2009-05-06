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
import math

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

##  Global Parameters of the simulation

SIZE    = 16.0  # size of the simulation box
NCHAINS = 5     # number of polymers
NBEADS  = 5     # number of particles in one polymer

class ParticleWriter(particles_PythonComputer):

    def __init__(self, _property, storage):
        particles_PythonComputer.__init__(self, storage)
        self.property = _property

    def each(self, pid):
        print("%d %s" % (pid, self.property[pid]))

import random

def randomWalk(step):

    rsq = 2.0
    while rsq > 1.0:
       dx = 2.0*random.random() - 1.0
       dy = 2.0*random.random() - 1.0
       dz = 2.0*random.random() - 1.0
       rsq = dx*dx + dy*dy + dz*dz

    r  = math.sqrt(rsq)

    return Real3D(step / r * dx, step / r * dy, step / r * dz)

particleStorage = Storage()
position = Real3DProperty(particleStorage)
velocity = Real3DProperty(particleStorage)
force    = Real3DProperty(particleStorage)

pbc = PBC()
pbc.set(SIZE)

# Bond List for bonded interactions

bondList = List(pbc, particleStorage, position)

for chainid in range(NCHAINS):
    pid1 = particleStorage.addParticle()
    pos1 = Real3D(random.random() * SIZE, random.random() * SIZE, random.random() * SIZE)
    position[pid1] = pos1
    velocity[pid1] = Real3D(0.0)
    force[pid1]    = Real3D(0.0)

    for beadid in range(NBEADS):

        # create a new position by random walk
        pos2 = pos1 + randomWalk(step=1.0)

        pid2 = particleStorage.addParticle()

        position[pid2] = pos2
        velocity[pid2] = Real3D(0.0)
        force[pid2]    = Real3D(0.0)

        bondList.addPair(pid1, pid2)

        # save pos2 and pid2 for the next step of this loop

        pos1 = pos2
        pid1 = pid2

#  Write all the particles of the storage

#  We could print the particles already here
#  particleStorage.foreach(ParticleWriter(position, particleStorage))

allSet = particles_All(particleStorage)
allSet.foreach(ParticleWriter(position, particleStorage))

# Pair List for nonbonded interactions

allPairs = All2(pbc, allSet, position)

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

