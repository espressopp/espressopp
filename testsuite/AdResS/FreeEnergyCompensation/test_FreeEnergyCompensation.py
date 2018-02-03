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
#

import sys
import time
import espressopp
import mpi4py.MPI as MPI

import unittest

class TestFreeEnergyCompensation(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        box = (10, 10, 10)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc=1.5,skin=system.skin)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc=1.5, skin=system.skin)
        system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)
        self.system = system

    def test_slab(self):
        # add some particles
        particle_list = [
            (1, 1, 0, espressopp.Real3D(5.5, 5.0, 5.0), 1.0, 0),
            (2, 1, 0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 0),
            (3, 1, 0, espressopp.Real3D(7.5, 5.0, 5.0), 1.0, 0),
            (4, 1, 0, espressopp.Real3D(8.5, 5.0, 5.0), 1.0, 0),
            (5, 1, 0, espressopp.Real3D(9.5, 5.0, 5.0), 1.0, 0),
            (6, 0, 0, espressopp.Real3D(5.5, 5.0, 5.0), 1.0, 1),
            (7, 0, 0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 1),
            (8, 0, 0, espressopp.Real3D(7.5, 5.0, 5.0), 1.0, 1),
            (9, 0, 0, espressopp.Real3D(8.5, 5.0, 5.0), 1.0, 1),
            (10, 0, 0, espressopp.Real3D(9.5, 5.0, 5.0), 1.0, 1),
        ]
        tuples = [(1,6),(2,7),(3,8),(4,9),(5,10)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()

        # generate a verlet list
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                dEx=2.0, dHy=1.0, adrCenter=[5.0, 5.0, 5.0], sphereAdr=False)

        # initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        # set up FEC
        fec = espressopp.integrator.FreeEnergyCompensation(self.system, center=[5.0, 5.0, 5.0])
        fec.addForce(itype=3,filename="table_fec.tab",type=1)
        integrator.addExtension(fec)

        # x coordinates of particles before integration
        before = [self.system.storage.getParticle(i).pos[0] for i in range(1,6)]

        # run ten steps and compute energy
        integrator.run(10)
        energy = fec.computeCompEnergy()

        # x coordinates of particles after integration
        after = [self.system.storage.getParticle(i).pos[0] for i in range(1,6)]

        # run checks (only one particle is in hybrid region and should feel the FEC. Also check that its FEC energy is correct)
        self.assertEqual(before[0], after[0])
        self.assertEqual(before[1], after[1])
        self.assertAlmostEqual(after[2], 7.598165, places=5)
        self.assertEqual(before[3], after[3])
        self.assertEqual(before[4], after[4])
        self.assertAlmostEqual(energy, 6.790157, places=5)

    def test_sphere(self):
        # add some particles
        particle_list = [
            (1, 1, 0, espressopp.Real3D(5.0, 5.5, 5.0), 1.0, 0),
            (2, 1, 0, espressopp.Real3D(5.0, 6.5, 5.0), 1.0, 0),
            (3, 1, 0, espressopp.Real3D(5.0, 7.5, 5.0), 1.0, 0),
            (4, 1, 0, espressopp.Real3D(5.0, 8.5, 5.0), 1.0, 0),
            (5, 1, 0, espressopp.Real3D(5.0, 9.5, 5.0), 1.0, 0),
            (6, 0, 0, espressopp.Real3D(5.0, 5.5, 5.0), 1.0, 1),
            (7, 0, 0, espressopp.Real3D(5.0, 6.5, 5.0), 1.0, 1),
            (8, 0, 0, espressopp.Real3D(5.0, 7.5, 5.0), 1.0, 1),
            (9, 0, 0, espressopp.Real3D(5.0, 8.5, 5.0), 1.0, 1),
            (10, 0, 0, espressopp.Real3D(5.0, 9.5, 5.0), 1.0, 1),
        ]
        tuples = [(1,6),(2,7),(3,8),(4,9),(5,10)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()

        # generate a verlet list
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                dEx=2.0, dHy=1.0, adrCenter=[5.0, 5.0, 5.0], sphereAdr=True)

        # initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        # set up FEC
        fec = espressopp.integrator.FreeEnergyCompensation(self.system, center=[5.0, 5.0, 5.0], sphereAdr=True)
        fec.addForce(itype=3,filename="table_fec.tab",type=1)
        integrator.addExtension(fec)

        # y coordinates of particles before integration
        before = [self.system.storage.getParticle(i).pos[1] for i in range(1,6)]

        # run ten steps
        integrator.run(10)
        energy = fec.computeCompEnergy()

        # y coordinates of particles after integration
        after = [self.system.storage.getParticle(i).pos[1] for i in range(1,6)]

        # run checks (as for test with slab-geometry, but check y-coordinates this time. Given the now spherical setup, particles should move as before but along the y-axis).
        self.assertEqual(before[0], after[0])
        self.assertEqual(before[1], after[1])
        self.assertAlmostEqual(after[2], 7.598165, places=5)
        self.assertEqual(before[3], after[3])
        self.assertEqual(before[4], after[4])
        self.assertAlmostEqual(energy, 6.790157, places=5)


if __name__ == '__main__':
    unittest.main()
