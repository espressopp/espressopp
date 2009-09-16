from espresso import *
from espresso import storage, bc, integrator, particles, potential, pairs, thermostat
import random

import logging

# logging.root.setLevel(logging.DEBUG)
# logging.getLogger("espresso.pmi").setLevel(logging.WARN)

size = 10.0
# set up the boundary conditions
pbc = bc.PeriodicBC(length=size)
# set up the decomposition scheme
storage = storage.SingleNode(pbc, 0)
#storage = espresso.storage.CellStorage(grid=(2, 2, 2), skin=0.1, **system)

# set up the properties
pos = storage.getPositionProperty()
vel = Real3DProperty(storage)
force = Real3DProperty(storage)
#force = storage.createProperty('Real3D')

# create the particles and the chains
bonds = pairs.List(bc=pbc, storage=storage)

# for chainid in range(100):
#     prevPid = storage.addParticle()
#     for beadid in range(63):
#         newPid = storage.addParticle()
#         pos[newPid] = pos[prevPid] \
#                       + Real3D(random.random()*size, random.random()*size, random.random()*size)
#         #+ espresso.esutil.randomWalk(step=1.0)
#         #bonds.add(particle, prevParticle)
#         prevPid = newPid

for beadid in range(30):
    storage.addParticle(beadid)
    pos[beadid] = Real3D(random.random()*size, random.random()*size, random.random()*size)
    # + espresso.esutil.randomWalk(step=1.0)
    vel[beadid] = 0.0
    force[beadid] = 0.0

pos[0] = (0.0, 0.0, 0.0)
pos[1] = (1.0, 0.0, 0.0)
bonds.addPair(0, 1)

# set up the integrator
integrator = integrator.VelocityVerlet(
    set=storage,
    velProperty=vel, forceProperty=force,
    timestep=0.01)

# create the LJ interaction
lj = potential.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.0)
if 0:
  # all pairs
  allpairs = pairs.All(set=storage)
  ljint = potential.Interaction(potential=lj, pairs=allpairs)
else:
  # Verlet lists
  vlist = pairs.VerletList(storage=storage, bc=pbc, radius=lj.cutoff, skin=0.3)
  ljint = potential.Interaction(potential=lj, pairs=vlist)
ljint.connect(integrator)

# create FENE interaction
fene = potential.FENE(K=1.0, r0=1.0, rMax=0.2)
feneint = potential.Interaction(potential=fene, pairs=bonds)
feneint.connect(integrator)

# set up the thermostat
therm = thermostat.Langevin(temperature=1.0, gamma=0.5)
therm.connect(integrator)

# # create the particles and the chains
# for chainid in range(100):
#     prevParticle = particles.addParticle()
#     for beadid in range(63):
#         newpos = lastParticle.pos + espresso.esutil.randomWalk(step=1.0)
#         particle = particles.addParticle(pos=newpos)
#         bonds.add(particle, prevParticle)
#         prevParticle = particle

# ke = analysis.KineticEnergy(storage, vel)

# calc the kinetic energy every 5 steps during integration

# integration

for sweeps in range(5):
    print('sweep ' + str(sweeps))
    print('  pos[0]=%8.3g %8.3g %8.3g' % tuple(pos[0]))
    print('  vel[0]=%8.3g %8.3g %8.3g' % tuple(vel[0]))
    print('  force[0]=%8.3g %8.3g %8.3g' % tuple(force[0]))

    print('  pos[1]=%8.3g %8.3g %8.3g' % tuple(pos[1]))
    print('  vel[1]=%8.3g %8.3g %8.3g' % tuple(vel[1]))
    print('  force[1]=%8.3g %8.3g %8.3g' % tuple(force[1]))
    integrator.integrate(steps=1)

#   print 'Kinetic energy = : ', ke.compute()
    print 'Pair energy = : ', ljint.totalEnergy()
    print 'Bond energy = : ', feneint.totalEnergy()
    print 'Pair virial = : ', ljint.totalVirial()
    print 'Bond virial = : ', feneint.totalVirial()

# analysis
