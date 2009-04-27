import espresso,
import espresso.bc
import espresso.decomposition
import espresso.thermostat
import espresso.integrator
import espresso.interaction
import espresso.pairs

system = {}
# set up the boundary conditions
system['bc'] = espresso.bc.PBC(length=10)
# set up the decomposition scheme
system['particles'] = espresso.decomposition.CellStorage(grid=(2, 2, 2), skin=0.1, **system)
# set up the thermostat
system['thermostat'] = espresso.thermostat.Langevin(temperatur=1.0, gamma=0.5, **system)
# set up the integrator
system['integrator'] = espresso.integrator.VelocityVerlet(timestep=0.001, **system)

# create the particles and the chains
feneint = espresso.interaction.FENE(K=1.0, r0=1.0, rMax=0.2)
bonds = espresso.pairs.List()
for chainid in range(100):
    prevParticle = particles.addParticle()
    for beadid in range(63):
        newpos = lastParticle.pos + espresso.esutil.randomWalk(step=1.0)
        particle = particles.addParticle(pos=newpos)
        bonds.add(particle, prevParticle)
        prevParticle = particle
system['integrator'].addInteraction(pairs=bonds, interaction=feneint)

# create the LJ interaction
ljint = espresso.interaction.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.0)
verletlists = espresso.pairs.VerletList(radius=ljint.cutoff, skin=0.3, **system)
system['integrator'].addInteraction(pairs=verletlists, interaction=ljint)

# integration
for sweeps in range(100):
    system['integrator'].integrate(steps=100)

# analysis
.
.
.
