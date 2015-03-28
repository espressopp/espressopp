#!/usr/bin/env python
"""
ESPresSo++
Test script for H5MD writer.
"""
import h5py
import os
import espressopp
import mpi4py.MPI as MPI
from espressopp.tools import decomp
from espressopp.io.DumpH5MD import *

box = (10.0, 10.0, 10.0)
rc = 2.5
skin = 0.3
system = espressopp.System()
system.rng = espressopp.esutil.RNG(12345)
system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin = skin

nodeGrid = decomp.nodeGrid(MPI.COMM_WORLD.size)
cellGrid = decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# Adds particles to the system and then decompose.
# id, type, mass, x, y, z
particles = [
    (1, 1, 1.0, espressopp.Real3D(1.0, 1.0, 1.0)),
    (2, 2, 1.0, espressopp.Real3D(2.0, 0.0, 2.0)),
    (3, 2, 1.0, espressopp.Real3D(3.0, 5.0, 1.0)),
    (4, 1, 1.0, espressopp.Real3D(6.0, 5.0, 4.0))
]
system.storage.addParticles(particles, 'id', 'type', 'mass', 'pos')
system.storage.decompose()

# integrator
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = 0.001

# h5md dump
h5md_group = 'atoms'
h5md_author = 'John Lenon'
h5md_email = 'john@lenon'
h5md_file = 'out.h5'
if os.path.exists(h5md_file):
    os.unlink(h5md_file)
with DumpH5MD(system, integrator, filename=h5md_file, h5md_group=h5md_group,
              author=h5md_author, email=h5md_email) as h5md:
    h5md.dump()

# Checks if the positions in h5md file are correct
h5 = h5py.File(h5md_file, 'r')
positions = list(h5['/particles/{}/position/value'.format(h5md_group)][0])
ids = list(h5['/particles/{}/id/value'.format(h5md_group)][0])
# Compare positions with particles
for p in range(len(positions)):
    assert list(positions[p]) == list(particles[p][3])

# Compare particle ids
for p in range(len(ids)):
    assert list(ids[p]) == [particles[p][0]]

# Clean up
os.unlink(h5md_file)

print '\nDone.'
