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


import sys
import time
import espressopp
import mpi4py.MPI as MPI

import unittest


class TestMinimizeEnergy(unittest.TestCase):
    def setUp(self):

        system = espressopp.System()
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, (10, 10, 10))
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        system.storage = espressopp.storage.DomainDecomposition(system)
        self.system = system


    def test_no_potential(self):
        particle_list = [
            (1, espressopp.Real3D(2.0, 2.0, 2.0), 1.0),
            (2, espressopp.Real3D(3.0, 2.0, 2.0), 1.0),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass')
        self.system.storage.decompose()
        minimize_energy = espressopp.integrator.MinimizeEnergy(self.system, 0.001, 0.0, 0.001)
        minimize_energy.run(10)
        self.assertEqual(minimize_energy.f_max, 0.0);

    def test_potential(self):
        particle_list = [
            (1, espressopp.Real3D(2.0, 2.0, 2.0), 1.0),
            (2, espressopp.Real3D(2.1, 2.0, 2.0), 1.0),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'pos', 'mass')
        self.system.storage.decompose()
        minimize_energy = espressopp.integrator.MinimizeEnergy(
            self.system, gamma=0.00001, ftol=1.0, max_displacement=0.001)

        vl = espressopp.VerletList(self.system, cutoff=2.5)
        lj = espressopp.interaction.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5)
        interaction = espressopp.interaction.VerletListLennardJones(vl)
        interaction.setPotential(type1=0, type2=0, potential=lj)
        self.system.addInteraction(interaction)

        energy_before = interaction.computeEnergy()

        minimize_energy.run(2000)
        self.assertLessEqual(minimize_energy.f_max, 1.0)
        self.assertLess(interaction.computeEnergy(), energy_before)


if __name__ == '__main__':
    unittest.main()
