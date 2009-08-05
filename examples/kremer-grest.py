from espresso import *
from espresso import storage, bc, integrator, particles, potential, pairs, thermostat
import random

size = 10.0
# set up the boundary conditions
pbc = bc.PeriodicBC(length=size)
# set up the decomposition scheme
storage = storage.SingleNode(0)
#storage = espresso.storage.CellStorage(grid=(2, 2, 2), skin=0.1, **system)

# set up the properties
pos = Real3DProperty(storage)
vel = Real3DProperty(storage)
force = Real3DProperty(storage)

# create the particles and the chains
feneint = potential.FENE(K=1.0, r0=1.0, rMax=0.2)
#bonds = pairs.List()


# for chainid in range(100):
#     prevPid = storage.addParticle()
#     for beadid in range(63):
#         newPid = storage.addParticle()
#         pos[newPid] = pos[prevPid] \
#                       + Real3D(random.random()*size, random.random()*size, random.random()*size)
#         #+ espresso.esutil.randomWalk(step=1.0)
#         #bonds.add(particle, prevParticle)
#         prevPid = newPid

for beadid in range(10):
    newPid = storage.addParticle()
    pos[newPid] = Real3D(random.random()*size, random.random()*size, random.random()*size)
        #+ espresso.esutil.randomWalk(step=1.0)
        #bonds.add(particle, prevParticle)
    vel[newPid] = 0.0
    force[newPid] = 0.0
    prevPid = newPid

# set up the integrator
integrator = integrator.VelocityVerlet(
    set=storage,
    posProperty=pos, velProperty=vel, forceProperty=force,
    timestep=0.01)

# create the LJ interaction
lj = potential.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.0)
allpairs = pairs.All(set=storage, bc=pbc, posProperty=pos)
ljint = potential.Interaction(potential=lj, pairs=allpairs)
ljint.connect(integrator)

# set up the thermostat
therm = thermostat.Langevin(temperature=1.0, gamma=0.5)
therm.connect(integrator)

# # create the particles and the chains
# feneint = espresso.potential.FENE(K=1.0, r0=1.0, rMax=0.2)
# bonds = espresso.pairs.List()
# for chainid in range(100):
#     prevParticle = particles.addParticle()
#     for beadid in range(63):
#         newpos = lastParticle.pos + espresso.esutil.randomWalk(step=1.0)
#         particle = particles.addParticle(pos=newpos)
#         bonds.add(particle, prevParticle)
#         prevParticle = particle
# system['integrator'].addInteraction(pairs=bonds, interaction=feneint)

#verletlists = espresso.pairs.VerletList(radius=ljint.cutoff, skin=0.3, **system)

# integration
for sweeps in range(5):
    print('sweep ' + str(sweeps))
    print('  pos[0]=%8.3g %8.3g %8.3g' % tuple(pos[0]))
    print('  vel[0]=%8.3g %8.3g %8.3g' % tuple(vel[0]))
    print('  force[0]=%8.3g %8.3g %8.3g' % tuple(force[0]))

    print('  pos[1]=%8.3g %8.3g %8.3g' % tuple(pos[1]))
    print('  vel[1]=%8.3g %8.3g %8.3g' % tuple(vel[1]))
    print('  force[1]=%8.3g %8.3g %8.3g' % tuple(force[1]))
    integrator.integrate(steps=1000)

# analysis
