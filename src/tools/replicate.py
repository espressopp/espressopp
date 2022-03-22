#  Copyright (C) 2012,2013,2016
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  Copyright (C) 2022
#      Max-Planck-Institute for Polymer Research & JGU Mainz
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

from espressopp import pmi
from espressopp import Real3D

class ReplicateParallelLocal:
    """ Usage:

        replicate   = (2,2,2)
        bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.lammps.read('polymer_melt.lammps')
        rp = espressopp.tools.ReplicateParallel()
        num_particles, Lx, Ly, Lz = rp.replicate(bonds, angles, x, y, z, Lx, Ly, Lz, *replicate)
        box = (Lx, Ly, Lz)
        system, integrator = espressopp.standard_system.Default(box=box,...)

        props = ['type', 'mass']
        num_particles_seed = len(x)
        seed_particles = []
        for i in range(num_particles_seed):
            part = [0, 1.0]
            seed_particles.append(part)

        rp.addParticles(system.storage, 1, seed_particles, *props)
        system.storage.decompose()

        fpl = espressopp.FixedPairList(system.storage)
        rp.addBonds(fpl)

        ftl = espressopp.FixedTripleList(system.storage)
        rp.addTriples(ftl)
    """

    def __init__(self):
        pass

    def replicateWorker(self, bonds, angles, x, y, z, Lx, Ly, Lz, xdim=1, ydim=1, zdim=1):
        self.bonds = bonds
        self.angles = angles
        self.x = x
        self.y = y
        self.z = z
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.xdim = xdim
        self.ydim = ydim
        self.zdim = zdim

    def addParticlesWorker(self, storage, start_pid, seed_particles, properties):
        assert len(seed_particles) == len(self.x), "Length mismatch in seed_particles"

        properties = ['id', 'pos', *properties]
        pid = start_pid
        for i in range(self.xdim):
            for j in range(self.ydim):
                for k in range(self.zdim):
                    new_particles = []
                    for x_, y_, z_, seed in zip(self.x, self.y, self.z, seed_particles):
                        pos = Real3D((x_ + i * self.Lx), (y_ + j * self.Ly), (z_ + k * self.Lz))
                        new_particles.append([pid, pos, *seed])
                        pid += 1
                    storage.addParticles(new_particles, *properties)

    def addBonds(self, fpl):
        ct = 0
        num_particles_original = len(self.x)
        fpl.addBonds(self.bonds)
        for i in range(self.xdim):
            for j in range(self.ydim):
                for k in range(self.zdim):
                    if(i + j + k != 0):
                        bonds_replicated = []
                        ct = ct + 1
                        for p1, p2 in self.bonds:
                            bonds_replicated.append((p1 + ct * num_particles_original, \
                                                    p2 + ct * num_particles_original))
                        fpl.addBonds(bonds_replicated)

    def addTriples(self, ftl):
        ct = 0
        num_particles_original = len(self.x)
        ftl.addTriples(self.angles)
        for i in range(self.xdim):
            for j in range(self.ydim):
                for k in range(self.zdim):
                    if(i + j + k != 0):
                        angles_replicated = []
                        for p1, p2, p3 in self.angles:
                            angles_replicated.append((p1 + ct * num_particles_original, \
                                                    p2 + ct * num_particles_original, \
                                                    p3 + ct * num_particles_original))
                        ftl.addTriples(angles_replicated)

if pmi.isController:
    class ReplicateParallel(metaclass=pmi.Proxy):
        pmiproxydefs = dict(
          cls = 'espressopp.tools.ReplicateParallelLocal',
          pmicall = ['replicateWorker', 'addParticlesWorker', 'addBonds', 'addTriples']
        )

        def replicate(self, bonds, angles, x, y, z, Lx, Ly, Lz, xdim=1, ydim=1, zdim=1):
            self.replicateWorker(bonds, angles, x, y, z, Lx, Ly, Lz, xdim, ydim, zdim)
            num_particles = len(x) * xdim * ydim * zdim
            Lx = Lx * xdim
            Ly = Ly * ydim
            Lz = Lz * zdim
            return num_particles, Lx, Ly, Lz

        def addParticles(self, storage, start_pid, seed_particles, *properties):
            if 'pos' in properties:
                raise ValueError("properties should not include 'pos'")
            if 'id' in properties:
                raise ValueError("properties should not include 'id'")
            self.addParticlesWorker(storage, start_pid, seed_particles, properties)

def replicate (bonds, angles, x, y, z, Lx, Ly, Lz, xdim=1, ydim=1, zdim=1):
    """
    Replicates configuration in each dimension.

    This may be used to increase the size of an equilibrated melt by a factor of 8 or more.

    Presently this routine works only for semiflexible polymers. A general
    class should be written to deal with files containing coordinates
    and topology data.

    xdim = ydim = zdim = 1 returns the original system not replicated.
    xdim = ydim = zdim = 2 returns the original system replicated to 8x.
    xdim = ydim = zdim = 3 returns the original system replicated to 27x.
    xdim = ydim = 1, zdim = 2 returns the original system replicated in the z-direction.
    """

    # replicate the particles
    x_replicated = x[:]
    y_replicated = y[:]
    z_replicated = z[:]
    for i in range(xdim):
        for j in range(ydim):
            for k in range(zdim):
                if(i + j + k != 0):
                    for x_, y_, z_ in zip(x, y, z):
                        x_replicated.append(x_ + i * Lx)
                        y_replicated.append(y_ + j * Ly)
                        z_replicated.append(z_ + k * Lz)

    # replicate the bonds and angles
    ct = 0
    num_particles_original = len(x)
    bonds_replicated = bonds[:]
    angles_replicated = angles[:]
    for i in range(xdim):
        for j in range(ydim):
            for k in range(zdim):
                if(i + j + k != 0):
                    ct = ct + 1
                    for p1, p2 in bonds:
                        bonds_replicated.append((p1 + ct * num_particles_original, \
                                                 p2 + ct * num_particles_original))
                    for p1, p2, p3 in angles:
                        angles_replicated.append((p1 + ct * num_particles_original, \
                                                  p2 + ct * num_particles_original, \
                                                  p3 + ct * num_particles_original))

    # modify the box size
    Lx = xdim * Lx
    Ly = ydim * Ly
    Lz = zdim * Lz

    return bonds_replicated, angles_replicated, x_replicated, y_replicated, z_replicated, Lx, Ly, Lz
