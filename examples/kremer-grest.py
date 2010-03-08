from espresso import *
from espresso import storage, bc, integrator, esutil
import random

# constants
size=10.0

# create objects
system = System()

# PROBLEM: The order is important, other order causes failure
system.rng = esutil.RNG()
system.bc = bc.OrthorhombicBC(system)
system.storage = storage.DomainDecomposition(system)

