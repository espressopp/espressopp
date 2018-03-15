#  Copyright (C) 2018
#      Max Planck Institute for Polymer Research
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

import sys
import time
import math
import espressopp
import mpi4py.MPI as MPI

import unittest

class TestCaseCapForce(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        rng = espressopp.esutil.RNG()
        rng.seed(1)
        system.rng = rng
        L = 10
        box = (L, L, L)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        self.system = system
        self.L = L
        self.box = box

    def test_cap_force(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(self.box, nodeGrid, 1.5, 0.3)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        particle_list = [
            (1, 0, 0, espressopp.Real3D(4.95, 5.0, 5.0), 1.0, 0, 1.),
            (2, 0, 0, espressopp.Real3D(5.05, 5.0, 5.0), 1.0, 0, 1.),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat', 'radius')
        self.system.storage.decompose()

        # integrator
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.005

        # Lennard-Jones with Verlet list
        rc_lj   = pow(2.0, 1.0/6.0)
        vl      = espressopp.VerletList(self.system, cutoff = rc_lj)
        potLJ   = espressopp.interaction.LennardJones(epsilon=1., sigma=1., cutoff=rc_lj, shift=0)
        interLJ = espressopp.interaction.VerletListLennardJones(vl)
        interLJ.setPotential(type1=0, type2=0, potential=potLJ)
        self.system.addInteraction(interLJ)

        # create a CapForce instance
        capforce       = espressopp.integrator.CapForce(self.system, 1.0)
        integrator.addExtension(capforce)

        # run 1 step
        integrator.run(1)

        particle1 = self.system.storage.getParticle(1)
        particle2 = self.system.storage.getParticle(2)
        print particle1.f, particle2.f

        # run checks
        self.assertTrue(math.fabs(particle1.f[0]) == 1.0, "The force of particle 1 is not capped.")
        self.assertTrue(math.fabs(particle1.f[0]) == 1.0, "The force of particle 2 is not capped.")

    def test_cap_force_array(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(self.box, nodeGrid, 1.5, 0.3)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        particle_list = [
            (1, 0, 0, espressopp.Real3D(4.95, 5.0, 5.0), 1.0, 0, 1.),
            (2, 0, 0, espressopp.Real3D(5.05, 5.0, 5.0), 1.0, 0, 1.),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat', 'radius')
        self.system.storage.decompose()

        # integrator
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.005

        # Lennard-Jones with Verlet list
        rc_lj   = pow(2.0, 1.0/6.0)
        vl      = espressopp.VerletList(self.system, cutoff = rc_lj)
        potLJ   = espressopp.interaction.LennardJones(epsilon=1., sigma=1., cutoff=rc_lj, shift=0)
        interLJ = espressopp.interaction.VerletListLennardJones(vl)
        interLJ.setPotential(type1=0, type2=0, potential=potLJ)
        self.system.addInteraction(interLJ)

        # create a CapForce instance
        capforce       = espressopp.integrator.CapForce(self.system, espressopp.Real3D(1.0, 1.0, 1.0))
        integrator.addExtension(capforce)

        # run 1 step
        integrator.run(1)

        particle1 = self.system.storage.getParticle(1)
        particle2 = self.system.storage.getParticle(2)
        print particle1.f, particle2.f

        # run checks
        self.assertTrue(math.fabs(particle1.f[0]) == 1.0, "The force of particle 1 is not capped.")
        self.assertTrue(math.fabs(particle1.f[0]) == 1.0, "The force of particle 2 is not capped.")

    def test_cap_force_group(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(self.box, nodeGrid, 1.5, 0.3)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        particle_list = [
            (1, 0, 0, espressopp.Real3D(4.95, 5.0, 5.0), 1.0, 0, 1.),
            (2, 0, 0, espressopp.Real3D(5.05, 5.0, 5.0), 1.0, 0, 1.),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat', 'radius')
        self.system.storage.decompose()

        # integrator
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.005

        # Lennard-Jones with Verlet list
        rc_lj   = pow(2.0, 1.0/6.0)
        vl      = espressopp.VerletList(self.system, cutoff = rc_lj)
        potLJ   = espressopp.interaction.LennardJones(epsilon=1., sigma=1., cutoff=rc_lj, shift=0)
        interLJ = espressopp.interaction.VerletListLennardJones(vl)
        interLJ.setPotential(type1=0, type2=0, potential=potLJ)
        self.system.addInteraction(interLJ)

        # create a ParticleGroup instance
        particle_group = espressopp.ParticleGroup(self.system.storage)
        particle_group.add(1)

        # create a CapForce instance
        capforce       = espressopp.integrator.CapForce(self.system, 1.0, particle_group)
        integrator.addExtension(capforce)

        # run 1 step
        integrator.run(1)

        particle1 = self.system.storage.getParticle(1)
        particle2 = self.system.storage.getParticle(2)
        print particle1.f, particle2.f

        # run checks
        self.assertTrue(math.fabs(particle1.f[0]) == 1.0, "The force of particle 1 is not capped.")
        self.assertTrue(math.fabs(particle2.f[0]) > 1.0, "The force of particle 2 is capped.")

    def test_cap_force_array_group(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(self.box, nodeGrid, 1.5, 0.3)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        particle_list = [
            (1, 0, 0, espressopp.Real3D(4.95, 5.0, 5.0), 1.0, 0, 1.),
            (2, 0, 0, espressopp.Real3D(5.05, 5.0, 5.0), 1.0, 0, 1.),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat', 'radius')
        self.system.storage.decompose()

        # integrator
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.005

        # Lennard-Jones with Verlet list
        rc_lj   = pow(2.0, 1.0/6.0)
        vl      = espressopp.VerletList(self.system, cutoff = rc_lj)
        potLJ   = espressopp.interaction.LennardJones(epsilon=1., sigma=1., cutoff=rc_lj, shift=0)
        interLJ = espressopp.interaction.VerletListLennardJones(vl)
        interLJ.setPotential(type1=0, type2=0, potential=potLJ)
        self.system.addInteraction(interLJ)

        # create a ParticleGroup instance
        particle_group = espressopp.ParticleGroup(self.system.storage)
        particle_group.add(1)

        # create a CapForce instance
        capforce       = espressopp.integrator.CapForce(self.system, espressopp.Real3D(1.0, 1.0, 1.0), particle_group)
        integrator.addExtension(capforce)

        # run 1 step
        integrator.run(1)

        particle1 = self.system.storage.getParticle(1)
        particle2 = self.system.storage.getParticle(2)
        print particle1.f, particle2.f

        # run checks
        self.assertTrue(math.fabs(particle1.f[0]) == 1.0, "The force of particle 1 is not capped.")
        self.assertTrue(math.fabs(particle2.f[0]) > 1.0, "The force of particle 2 is capped.")

if __name__ == '__main__':
    unittest.main()
