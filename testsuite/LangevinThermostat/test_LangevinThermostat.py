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

class TestLangevinThermostat(unittest.TestCase):
    global box
    box = (10, 10, 10)		
    def setUp(self):
        # set up system
        system = espressopp.System()
        rng = espressopp.esutil.RNG()
        rng.seed(1)
        system.rng = rng
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        self.system = system

    def test_normal(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc=1.5,skin=0.3)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc=1.5, skin=0.3)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        particle_list = [
            (1, 1, 0, espressopp.Real3D(5.5, 5.0, 5.0), 1.0, 0),
            (2, 1, 0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 0),
            (3, 1, 0, espressopp.Real3D(7.5, 5.0, 5.0), 1.0, 0),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat')
        self.system.storage.decompose()

        # generate a verlet list
        vl = espressopp.VerletList(self.system, cutoff=1.5)

        # integrator
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01

        # Langevin Thermostat
        langevin = espressopp.integrator.LangevinThermostat(self.system)
        langevin.gamma = 1.0
        langevin.temperature = 1.0
        langevin.adress = False
        langevin.addExclusions([1])
        integrator.addExtension(langevin)

        # coordinates of particles before integration
        before =[self.system.storage.getParticle(i).pos[j] for i in range(1,4) for j in range(3)]

        # run ten steps
        integrator.run(10)

        # coordinates of particles after integration
        after = [self.system.storage.getParticle(i).pos[j] for i in range(1,4) for j in range(3)]

        # run checks (first particle excluded, hence it should not move. The other should have moved, however, as they feel the thermostat)
        self.assertEqual(before[0], after[0])
        self.assertEqual(before[1], after[1])
        self.assertEqual(before[2], after[2])
        self.assertNotEqual(before[3], after[3])
        self.assertNotEqual(before[4], after[4])
        self.assertNotEqual(before[5], after[5])
        self.assertNotEqual(before[6], after[6])
        self.assertNotEqual(before[7], after[7])
        self.assertNotEqual(before[8], after[8])

    def test_AdResS(self):
        # set up AdResS domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc=1.5,skin=0.3)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc=1.5, skin=0.3)
        self.system.storage = espressopp.storage.DomainDecompositionAdress(self.system, nodeGrid, cellGrid)

        # add some particles (atomistic and coarse-grained particles now)
        particle_list = [
            (1, 1, 0, espressopp.Real3D(5.5, 5.0, 5.0), 1.0, 0),
            (2, 1, 0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 0),
            (3, 1, 0, espressopp.Real3D(7.5, 5.0, 5.0), 1.0, 0),
            (4, 0, 0, espressopp.Real3D(5.5, 5.0, 5.0), 1.0, 1),
            (5, 0, 0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 1),
            (6, 0, 0, espressopp.Real3D(7.5, 5.0, 5.0), 1.0, 1),
        ]
        tuples = [(1,4),(2,5),(3,6)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()

        # generate a verlet list
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                dEx=1.0, dHy=1.0, adrCenter=[5.0, 5.0, 5.0], sphereAdr=False)

        # integrator including AdResS
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system, vl, ftpl)
        integrator.addExtension(adress)

        # Langevin Thermostat
        langevin = espressopp.integrator.LangevinThermostat(self.system)
        langevin.gamma = 1.0
        langevin.temperature = 1.0
        langevin.adress = True
        langevin.addExclusions([4])
        integrator.addExtension(langevin)

        # coordinates of particles before integration
        before = [self.system.storage.getParticle(i).pos[j] for i in range(1,4) for j in range(3)]

        # run ten steps
        integrator.run(10)

        # coordinates of particles after integration
        after = [self.system.storage.getParticle(i).pos[j] for i in range(1,4) for j in range(3)]

        # run checks (same as test before)
        self.assertEqual(before[0], after[0])
        self.assertEqual(before[1], after[1])
        self.assertEqual(before[2], after[2])
        self.assertNotEqual(before[3], after[3])
        self.assertNotEqual(before[4], after[4])
        self.assertNotEqual(before[5], after[5])
        self.assertNotEqual(before[6], after[6])
        self.assertNotEqual(before[7], after[7])
        self.assertNotEqual(before[8], after[8])


if __name__ == '__main__':
    unittest.main()
