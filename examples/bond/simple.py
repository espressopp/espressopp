import espresso
import espresso.bc
import espresso.pairs
import random

from espresso import Real3D
from espresso import Real3DProperty
from espresso import storage
from espresso import integrator
from espresso import potential

import logging 
logging.root.setLevel(logging.DEBUG)
logging.getLogger("espresso.pmi").setLevel(logging.WARN)

size = 10.0
pbc=espresso.bc.PeriodicBC(length=10.0)
storage = storage.SingleNode(pbc, 0)

# set up the properties
pos   = Real3DProperty(storage)
vel   = Real3DProperty(storage)
force = Real3DProperty(storage)

for pid in range(10):
    storage.addParticle(pid)
    pos[pid] = Real3D(random.random()*size, random.random()*size, random.random()*size)
    vel[pid] = Real3D(random.random(), random.random(), random.random())
    # vel[pid] = Real3D(0.0)

allpairs=espresso.pairs.All(set=storage)

# default: shift=auto, offset=0.0
ljpot=espresso.potential.LennardJones(sigma=1, epsilon=1, cutoff=2.5)

# Just compute the forces
# allpairs.computeForces(set='force', interaction=ljint)
# force=particles.getProperty(name='force')
# print(force[particles[17]])
# print(particles[17].get(name='force'))

# Now do  simulation steps

vvintegrator=integrator.VelocityVerlet(set = storage,
    velProperty = vel, forceProperty = force, timestep = 1.0)

ljint = espresso.potential.Interaction(ljpot, allpairs)

ljint.connect(vvintegrator)

vvintegrator.integrate(steps=1)

print(force[7])
