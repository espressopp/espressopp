#!/usr/bin/env python3

import espressopp
from espressopp import Real3D, Int3D
from espressopp.tools import decomp, lattice, velocities

from mpi4py import MPI
import time
import numpy as np

def generate_particles(particles_per_direction):
    num_particles = particles_per_direction**3
    x, y, z, Lx, Ly, Lz = lattice.createCubic(num_particles, rho=0.8442, perfect=False)
    vx, vy, vz = velocities.gaussian(T=0.6, N=num_particles, zero_momentum=True)
    return x, y, z, Lx, Ly, Lz, vx, vy, vz

def generate_system(add_particles_array):
    rc = 2.5
    skin = 0.3
    timestep = 0.005
    temperature = 1.0
    comm = MPI.COMM_WORLD

    particles_per_direction = 64
    x, y, z, Lx, Ly, Lz, vx, vy, vz = generate_particles(particles_per_direction)

    num_particles = len(x)
    density = num_particles / (Lx * Ly * Lz)
    size = (Lx, Ly, Lz)

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG()
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
    system.skin = skin
    nodeGrid = decomp.nodeGrid(comm.size, size, rc, skin)
    cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

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

    return num_particles, tprep, tadd

if __name__=="__main__":
    num_particles1, tprep1, tadd1 = generate_system(False)
    num_particles2, tprep2, tadd2 = generate_system(True)

    print("\n")
    print("Using system.storage.addParticles(...)")
    print("    prepared {} particles in {:.2f} seconds".format(num_particles1, tprep1))
    print("    added    {} particles in {:.2f} seconds".format(num_particles1, tadd1))
    print()
    print("Using system.storage.addParticlesArray(...)")
    print("    prepared {} particles in {:.2f} seconds, speedup: {:.2f}".format(num_particles2, tprep2, tprep1/tprep2))
    print("    added    {} particles in {:.2f} seconds, speedup: {:.2f}".format(num_particles2, tadd2, tadd1/tadd2))
