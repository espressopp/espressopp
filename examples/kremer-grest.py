import espresso,
import espresso.bc
import espresso.decomposition
import espresso.thermostat
import espresso.integrator
import espresso.interaction
import espresso.pairs

# set up the boundary conditions
pbc = espresso.bc.PBC(length=10)
# set up the decomposition scheme
particles = espresso.decomposition.CellStorage(bc=pbc, grid=(2, 2, 2), skin=0.1)
# set up the thermostat
lvThermostat = espresso.thermostat.Langevin(temperatur=1.0, gamma=0.5)
# set up the integrator
vvintegrator = espresso.integrator.VelocityVerlet(
    timestep=0.001, thermostat=lvThermostat)

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
vvintegrator.addInteraction(pairs=bonds, interaction=feneint)

# create the LJ interaction
ljint = espresso.interaction.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.0)
verletlists = espresso.pairs.VerletList(radius=ljint.cutoff, skin=0.3)
vvintegrator.addInteraction(pairs=verletlists, interaction=ljint)

# integration
for sweeps in range(100):
    vvintegrator.integrate(steps=100)

# analysis
.
.
.
