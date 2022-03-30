import unittest, time
import numpy as np
from mpi4py import MPI
import espressopp
from espressopp import Real3D, Int3D
from espressopp.tools import decomp, lattice, velocities

def generate_particles(particles_per_direction):
    num_particles = particles_per_direction**3
    x, y, z, Lx, Ly, Lz = lattice.createCubic(
        num_particles, rho=0.8442, perfect=False)
    vx, vy, vz = velocities.gaussian(
        T=0.6, N=num_particles, zero_momentum=True)
    return x, y, z, Lx, Ly, Lz, vx, vy, vz

particles_per_direction = 10
partices_orig = generate_particles(particles_per_direction)

def generate_system(Lx, Ly, Lz, num_particles):
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

def add_particles(add_particles_array, repl=None):

    x, y, z, Lx, Ly, Lz, vx, vy, vz = partices_orig
    num_particles_orig = len(x)

    if repl is not None:
        box_orig = (Lx, Ly, Lz)
        num_particles, Lx, Ly, Lz = espressopp.storage.preReplicate(num_particles_orig, *box_orig, *repl)
    else:
        num_particles = num_particles_orig

    system = generate_system(Lx, Ly, Lz, num_particles)

    if add_particles_array:

        tstart = time.time()
        props = ['id', 'type', 'mass', 'posx',
                 'posy', 'posz', 'vx', 'vy', 'vz']
        ids = np.arange(1, num_particles_orig+1)
        types = np.zeros(num_particles_orig)
        mass = np.ones(num_particles_orig)
        new_particles = np.stack(
            (ids, types, mass, x, y, z, vx, vy, vz), axis=-1)
        tprep = time.time()-tstart

        if repl is None:
            tstart = time.time()
            system.storage.addParticlesArray(new_particles, *props)
            tadd = time.time()-tstart
        else:
            tstart = time.time()
            system.storage.addParticlesArrayReplicate(
                new_particles, *props, size_orig=(Lx, Ly, Lz), repl=repl)
            tadd = time.time()-tstart

    else:

        if repl is not None:
            __b, __a, x, y, z, Lx, Ly, Lz = espressopp.tools.replicate([], [], x, y, z, Lx, Ly, Lz,
                                                                       xdim=repl[0], ydim=repl[1], zdim=repl[2])
            tot_repl = repl[0]*repl[1]*repl[2]
            vx = np.tile(vx, tot_repl)
            vy = np.tile(vy, tot_repl)
            vz = np.tile(vz, tot_repl)

        tstart = time.time()
        props = ['id', 'type', 'mass', 'pos', 'v']
        new_particles = []
        for i in range(num_particles):
            new_particles.append(
                [i + 1, 0, 1.0, Real3D(x[i], y[i], z[i]), Real3D(vx[i], vy[i], vz[i])])
        tprep = time.time()-tstart

        tstart = time.time()
        system.storage.addParticles(new_particles, *props)
        tadd = time.time()-tstart

    system.storage.decompose()

    return system

def get_particles_config(system):

    tstart = time.time()
    configurations = espressopp.analysis.Configurations(system, pos=True, vel=True)
    configurations.gather()
    conf = configurations[0]
    tconfig = time.time()-tstart
    print(f"tconfig = {tconfig} s")

    ids = conf.getIds()
    pos = [conf.getCoordinates(pid) for pid in ids]
    vel = [conf.getVelocities(pid) for pid in ids]

    return pos, vel

def get_particles_array(system):
    tstart = time.time()
    properties = ["id", "posx", "posy", "posz", "vx", "vy", "vz", "type", "mass"]
    data  = system.storage.getParticlesArray(*properties)
    tarray = time.time()-tstart
    print(f"tarray = {tarray} s")

    pid = data[0]
    pid_set = set(pid)
    assert len(pid)==len(pid_set)

    pid_sort = pid.argsort()

    pos = data[1:4].transpose()[pid_sort]
    vel = data[4:7].transpose()[pid_sort]

    return pos, vel

def add_missing_props():
    x, y, z, Lx, Ly, Lz, vx, vy, vz = partices_orig
    num_particles = len(x)
    system = generate_system(Lx, Ly, Lz, num_particles)
    props = ['id', 'type', 'mass', 'posx', 'posy', 'posz', 'vx', 'vy']
    ids = np.arange(1, num_particles+1)
    types = np.zeros(num_particles)
    mass = np.ones(num_particles)
    new_particles = np.stack((ids, types, mass, x, y, z, vx, vy, vz), axis=-1)

    system.storage.addParticlesArray(new_particles, *props)
