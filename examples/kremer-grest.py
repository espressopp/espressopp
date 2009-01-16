import espresso

# - setup simulation system
#    - set system geometry
#    - set integrator
#    - set particles
#    - set interactions
#      - set bonded interactions
#      - set nonbonded interactions
#    - set decomposition type
# - setup topology
# - integrate
# - analysis

#####
system = espresso.system.System()
# default: type=float, dim=1
velocity = system.addParticleProperty(name="velocity", dim=3)

## SET GEOMETRY
system.geometry = espresso.bc.PBC(length = 10)

## SET INTEGRATOR
lvthermostat = espresso.LangevinThermostat(temperature=1.0, gamma=0.5)
# default: useAsVelocity="velocity", useAsForce="force"
vvintegrator = espresso.VelocityVerletIntegrator(timestep=0.001)
# default: thermostat=none
vvintegrator.thermostat = lvthermostat
system.integrator = vvintegrator

## SET DECOMPOSITION
# default: grid, skin
domdec = espresso.DomainDecomposition(grid = 2 2 2, skin = 0.1)
system.decomposition = domdec

## SET CHAINS
# default: dim=2
chains = espresso.ParticleTuple(dim=2)
for chainid in range(100):
    # default: pos=random
    lastParticle = system.addParticle()
    for beadid in range(63):
        newpos = lastParticle.pos + randomWalk(step=1.0)
        particle = system.addParticle(pos=newpos)
        chains.addTuple(particle, lastParticle)
        lastParticle = particle

## SET INTERACTIONS
# default: FENE parameters
system.addInteraction(tuples=chains, interaction=espresso.FENEInteraction(...))

# default: particles=all, exclusions=none
verletlists = espresso.VerletListsTuples(skin=0.3)
# default: shift=auto, offset=0.0
ljint = espresso.LJPairInteraction(sigma=1.0, epsilon=1.0, cutoff=1.0)
system.addInteraction(tuples=verletlists, tupleInteraction=ljint)

## INTEGRATE
for sweeps in range(100):
    system.integrate(steps=100)
    # - write system to disc
    # - VMD coupling
    # - analysis

# example for selecting a partilcle group:
pg = system.selectParticles("within 5 of 3.0 3.0 3.0")
