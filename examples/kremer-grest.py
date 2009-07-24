from espresso import *
from espresso import bc, decomposition, integrator, particles, potential, pairs

# import espresso.thermostat
# import espresso.interaction
# import espresso.pairs

# set up the boundary conditions
pbc = bc.PeriodicBC(length=10)
# set up the decomposition scheme
decomposer = decomposition.SingleNode(0)
#particles = espresso.decomposition.CellStorage(grid=(2, 2, 2), skin=0.1, **system)

# set up the properties
pos = decomposer.createProperty('Real3D')
vel = decomposer.createProperty('Real3D')
force = decomposer.createProperty('Real3D')

# create the particles and the chains
feneint = potential.FENE(K=1.0, r0=1.0, rMax=0.2)
#bonds = pairs.List()
for chainid in range(100):
    prevPid = decomposer.addParticle()
    for beadid in range(63):
        newPid = decomposer.addParticle()
        pos[newPid] = pos[prevPid]
        #+ espresso.esutil.randomWalk(step=1.0)
        #bonds.add(particle, prevParticle)
        prevPid = newPid

# create the LJ interaction
lj = potential.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.0)
allpairs = pairs.All(set=decomposer, bc=pbc, posProperty=pos)
ljint = potential.Interaction(potential=lj, pairs=allpairs)

# set up the integrator
integrator = integrator.VelocityVerlet(
    set=decomposer,
    posProperty=pos, velProperty=vel, forceProperty=force,
    timestep=0.001)
# set up the thermostat
thermostat = thermostat.Langevin(temperature=1.0, gamma=0.5)

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
# for sweeps in range(100):
#     integrator.integrate(steps=100)

# analysis
