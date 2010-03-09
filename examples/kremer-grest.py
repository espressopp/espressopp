from espresso import *
from espresso import storage, bc, integrator, esutil, interaction
import random

# constants
size=10.0

# create objects
system = System()

# PROBLEM: The order is important, other order causes failure
system.rng = esutil.RNG()
system.bc = bc.OrthorhombicBC(system)
system.storage = storage.DomainDecomposition(system)

vl = VerletList(system, cutoff=1.0, skin=0.4)
ljint = interaction.VerletListLennardJones(vl)
ljfunc = interaction.LennardJonesFunction(
    epsilon=1.0, sigma=1.0, cutoff=2.0)
ljint.setParameters(type1=0, type2=1, ljfunc)

integrator = integrator.VelocityVerlet(system, timestep=0.001)
integrator.addInteraction(ljint)

integrator.run(steps=1000)
