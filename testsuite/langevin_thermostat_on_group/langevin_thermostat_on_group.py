#!/usr/bin/env python2
#
#  Copyright (C) 2013-2017(H)
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
# 
# -*- coding: utf-8 -*-


import espressopp
import mpi4py.MPI as MPI

import unittest

class TestCaseLangevinThermostatOnGroup(unittest.TestCase):
    def setUp(self):
        system = espressopp.System()
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, (10, 10, 10))
        system.rng = espressopp.esutil.RNG(54321)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        self.system = system

        system.storage = espressopp.storage.DomainDecomposition(system)

        particle_lists = [
            (1, espressopp.Real3D(0, 1, 2), espressopp.Real3D(0, 0, 0)),
            (2, espressopp.Real3D(0, 1, 2), espressopp.Real3D(0, 0, 0)),
            (3, espressopp.Real3D(0, 1, 2), espressopp.Real3D(0, 0, 0)),
            (4, espressopp.Real3D(0, 1, 2), espressopp.Real3D(0, 0, 0)),
            (5, espressopp.Real3D(0, 1, 2), espressopp.Real3D(0, 0, 0))
        ]
        self.thermo_group_pids = [1, 2, 3]
        self.non_thermo_group_pids = [4, 5]
        system.storage.addParticles(particle_lists, 'id', 'pos', 'v')

        self.thermo_group = espressopp.ParticleGroup(system.storage)

        for p in self.thermo_group_pids:
            self.thermo_group.add(p)

        self.integrator = espressopp.integrator.VelocityVerlet(system)
        self.integrator.dt = 0.001

    def test_thermalize_only_group(self):
        """Run thermostat only on particles in ParticleGroup."""
        langevin = espressopp.integrator.LangevinThermostatOnGroup(self.system, self.thermo_group)
        langevin.gamma = 1.0
        langevin.temperature = 1.0
        self.integrator.addExtension(langevin)

        self.integrator.run(1)

        # Compare the forces on particles, 4 and 5 == 0
        forces = [self.system.storage.getParticle(p).f for p in range(1, 6)]
        for p in self.thermo_group_pids:
            f = self.system.storage.getParticle(p).f
            self.assertNotEqual(f, espressopp.Real3D())

        for p in self.non_thermo_group_pids:
            f = self.system.storage.getParticle(p).f
            self.assertEqual(f, espressopp.Real3D())


if __name__ == '__main__':
    unittest.main()
