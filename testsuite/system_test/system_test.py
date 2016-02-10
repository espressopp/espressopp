#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import time
import espressopp
import mpi4py.MPI as MPI
from espressopp.tools import timers


system = espressopp.System()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, (10, 10, 10))
system.skin = 0.3
system.comm = MPI.COMM_WORLD
#nodeGrid = decomp.nodeGrid(comm.size)
#cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system)

vl = espressopp.VerletList(system, cutoff=1.4)
lj1 = espressopp.interaction.VerletListLennardJones(vl)
lj2 = espressopp.interaction.VerletListLennardJones(vl)
lj3 = espressopp.interaction.VerletListLennardJones(vl)
lj4 = espressopp.interaction.VerletListLennardJones(vl)

system.addInteraction(lj1, 'lj1')
system.addInteraction(lj2, 'lj2')
system.addInteraction(lj3, 'lj3')
system.addInteraction(lj4)

# Checks if all interactions are present.
assert set(['lj1', 'lj2', 'lj3']) == set(system.getAllInteractions().keys())
assert id(lj1) == id(system.getInteractionByName('lj1'))

assert id(lj1) == id(system.getInteraction(0))
assert id(lj2) == id(system.getInteraction(1))
assert id(lj3) == id(system.getInteraction(2))

# Checks remove by name.
assert system.getNumberOfInteractions() == 4
system.removeInteractionByName('lj2')
assert system.getNumberOfInteractions() == 3

# Checks if interactions are renumbered correctly.
assert id(lj1) == id(system.getInteractionByName('lj1'))
assert id(lj3) == id(system.getInteractionByName('lj3'))
