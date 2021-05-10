#!/usr/bin/env python3
#  Copyright (C) 2021
#      Max Planck Institute for Polymer Research & JGU Mainz
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest, time
import numpy as np
from mpi4py import MPI
import espressopp
from espressopp import Real3D, Int3D
from espressopp.tools import decomp, lattice, velocities

def generate_particles(particles_per_direction):
    num_particles = particles_per_direction**3
    x, y, z, Lx, Ly, Lz = lattice.createCubic(num_particles, rho=0.8442, perfect=False)
    vx, vy, vz = velocities.gaussian(T=0.6, N=num_particles, zero_momentum=True)
    return x, y, z, Lx, Ly, Lz, vx, vy, vz

particles_per_direction = 10
x, y, z, Lx, Ly, Lz, vx, vy, vz = generate_particles(particles_per_direction)
num_particles = len(x)

def generate_system():
    rc = 2.5
    skin = 0.3
    timestep = 0.005
    temperature = 1.0
    comm = MPI.COMM_WORLD
    density = num_particles / (Lx * Ly * Lz)
    size = (Lx, Ly, Lz)

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG()
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
    system.skin = skin
    nodeGrid = decomp.nodeGrid(comm.size, size, rc, skin)
    cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    return system

def add_particles(add_particles_array):

    system = generate_system()

    if add_particles_array:

        tstart = time.time()
        props = ['id', 'type', 'mass', 'posx', 'posy', 'posz', 'vx', 'vy', 'vz']
        ids      = np.arange(1,num_particles+1)
        types    = np.zeros(num_particles)
        mass     = np.ones(num_particles)
        new_particles = np.stack((ids, types, mass, x, y, z, vx, vy, vz), axis=-1)
        tprep = time.time()-tstart

        tstart = time.time()
        system.storage.addParticlesArray(new_particles, *props)
        tadd = time.time()-tstart

    else:

        tstart = time.time()
        props = ['id', 'type', 'mass', 'pos', 'v']
        new_particles = []
        for i in range(num_particles):
            new_particles.append([i + 1, 0, 1.0, Real3D(x[i], y[i], z[i]), Real3D(vx[i], vy[i], vz[i])])
        tprep = time.time()-tstart

        tstart = time.time()
        system.storage.addParticles(new_particles, *props)
        tadd = time.time()-tstart

    system.storage.decompose()

    configurations = espressopp.analysis.Configurations(system, pos=True, vel=True)
    configurations.gather()
    conf = configurations[0]

    ids = conf.getIds()
    pos = [conf.getCoordinates(pid) for pid in ids]
    vel = [conf.getVelocities(pid) for pid in ids]

    return pos, vel

def add_missing_props():
    system = generate_system()
    props = ['id', 'type', 'mass', 'posx', 'posy', 'posz', 'vx', 'vy']
    ids      = np.arange(1,num_particles+1)
    types    = np.zeros(num_particles)
    mass     = np.ones(num_particles)
    new_particles = np.stack((ids, types, mass, x, y, z, vx, vy, vz), axis=-1)

    system.storage.addParticlesArray(new_particles, *props)

class TestAddParticlesArray(unittest.TestCase):

    def compare(self, prop0, prop1):
        self.assertEqual(len(prop0), len(prop1))
        for i in range(len(prop0)):
            for j in range(3):
                self.assertEqual(prop0[i][j],prop1[i][j])

    def test1(self):
        pos0, vel0 = add_particles(False)
        pos1, vel1 = add_particles(True)

        self.compare(pos0, pos1)
        self.compare(vel0, vel1)

    def test2(self):
        with self.assertRaises(AssertionError):
            add_missing_props()

if __name__ == "__main__":
    unittest.main()
