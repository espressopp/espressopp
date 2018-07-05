#!/usr/bin/env python2
#
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
#
# -*- coding: utf-8 -*-

import espressopp
import mpi4py.MPI as MPI

import unittest

class TestMTSAdResS(unittest.TestCase):
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

        # add some particles
        particle_list = [
            (1, 1, 0, espressopp.Real3D(5.5, 5.0, 5.0), 1.0, 0),
            (2, 1, 0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 0),
            (3, 1, 0, espressopp.Real3D(7.5, 5.0, 5.0), 1.0, 0),
            (4, 1, 0, espressopp.Real3D(8.5, 5.0, 5.0), 1.0, 0),
            (5, 1, 0, espressopp.Real3D(9.5, 5.0, 5.0), 1.0, 0),
            (6, 0, 1.0, espressopp.Real3D(5.5, 5.0, 5.0), 1.0, 1),
            (7, 0, 1.0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 1),
            (8, 0, 1.0, espressopp.Real3D(7.5, 5.0, 5.0), 1.0, 1),
            (9, 0, 1.0, espressopp.Real3D(8.5, 5.0, 5.0), 1.0, 1),
            (10, 0, 1.0, espressopp.Real3D(9.5, 5.0, 5.0), 1.0, 1),
        ]
        tuples = [(1,6),(2,7),(3,8),(4,9),(5,10)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat')
        self.ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        self.ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(self.ftpl)
        self.system.storage.decompose()

        # set up a verlet list
        self.vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                dEx=2.0, dHy=1.0, adrCenter=[5.0, 5.0, 5.0], sphereAdr=False)

        # set up integrator
        self.integrator = espressopp.integrator.VelocityVerletRESPA(self.system)
        self.integrator.dt = 0.01
        self.integrator.multistep = 4
        adress = espressopp.integrator.Adress(self.system, self.vl, self.ftpl, multistep=4)
        self.integrator.addExtension(adress)

    def test_ATAT_CG_templates_HAdResS(self):
        # add atomistic interactions
        interNBat = espressopp.interaction.VerletListHadressATLenJonesReacFieldGen(self.vl, self.ftpl)
        potLJ  = espressopp.interaction.LennardJones(epsilon=0.650299305951, sigma=0.316549165245, shift='auto', cutoff=1.4)
        potQQ  = espressopp.interaction.ReactionFieldGeneralized(prefactor=138.935485, kappa=0.0, epsilon1=1.0, epsilon2=80.0, cutoff= 1.4, shift="auto")
        interNBat.setPotential1(type1=0, type2=0, potential=potLJ)
        interNBat.setPotential2(type1=0, type2=0, potential=potQQ)
        self.system.addInteraction(interNBat)

        # add coarse-grained interaction
        interNBcg = espressopp.interaction.VerletListHadressCGHarmonic(self.vl, self.ftpl)
        potCG = espressopp.interaction.Harmonic(K=500.0, r0=1.4, cutoff=1.4)
        interNBcg.setPotential(type1=1, type2=1, potential=potCG)
        self.system.addInteraction(interNBcg)

        # decompose properly
        espressopp.tools.AdressDecomp(self.system, self.integrator)

        # coordinates and non-bonded energy of particles before integration
        before = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_before = interNBcg.computeEnergy() + interNBat.computeEnergy()

        # run ten steps and compute energy
        self.integrator.run(10)

        # coordinates and non-bonded energy of particles after integration
        after = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_after = interNBcg.computeEnergy() + interNBat.computeEnergy()

        # run checks
        self.assertAlmostEqual(after[0], 3.513552, places=5)
        self.assertEqual(before[1], after[1])
        self.assertEqual(before[2], after[2])
        self.assertAlmostEqual(after[3], 4.260359, places=5)
        self.assertEqual(before[4], after[4])
        self.assertEqual(before[5], after[5])
        self.assertAlmostEqual(after[6], 7.384342, places=5)
        self.assertEqual(before[7], after[7])
        self.assertEqual(before[8], after[8])
        self.assertAlmostEqual(after[9], 9.690175, places=5)
        self.assertEqual(before[10], after[10])
        self.assertEqual(before[11], after[11])
        self.assertAlmostEqual(after[12], 1.239614, places=5)
        self.assertEqual(before[13], after[13])
        self.assertEqual(before[14], after[14])
        self.assertAlmostEqual(energy_before, 191.685729, places=5)
        self.assertAlmostEqual(energy_after, 51.946263, places=5)

    def test_AT_CG_templates_HAdResS(self):
        # add atomistic interaction
        interNBat = espressopp.interaction.VerletListHadressATLennardJones(self.vl, self.ftpl)
        potLJ  = espressopp.interaction.LennardJones(epsilon=0.650299305951, sigma=0.316549165245, shift='auto', cutoff=1.4)
        interNBat.setPotential(type1=0, type2=0, potential=potLJ)
        self.system.addInteraction(interNBat)

        # add coarse-grained interaction
        interNBcg = espressopp.interaction.VerletListHadressCGHarmonic(self.vl, self.ftpl)
        potCG = espressopp.interaction.Harmonic(K=50.0, r0=1.4, cutoff=1.4)
        interNBcg.setPotential(type1=1, type2=1, potential=potCG)
        self.system.addInteraction(interNBcg)

        # decompose properly
        espressopp.tools.AdressDecomp(self.system, self.integrator)

        # coordinates and non-bonded energy of particles before integration
        before = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_before = interNBcg.computeEnergy() + interNBat.computeEnergy()

        # run ten steps and compute energy
        self.integrator.run(10)

        # coordinates and non-bonded energy of particles after integration
        after = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_after = interNBcg.computeEnergy() + interNBat.computeEnergy()

        # run checks
        self.assertAlmostEqual(after[0], 5.507295, places=5)
        self.assertEqual(before[1], after[1])
        self.assertEqual(before[2], after[2])
        self.assertAlmostEqual(after[3], 5.825483, places=5)
        self.assertEqual(before[4], after[4])
        self.assertEqual(before[5], after[5])
        self.assertAlmostEqual(after[6], 6.408025, places=5)
        self.assertEqual(before[7], after[7])
        self.assertEqual(before[8], after[8])
        self.assertAlmostEqual(after[9], 8.368686, places=5)
        self.assertEqual(before[10], after[10])
        self.assertEqual(before[11], after[11])
        self.assertAlmostEqual(after[12], 0.763036, places=5)
        self.assertEqual(before[13], after[13])
        self.assertEqual(before[14], after[14])
        self.assertAlmostEqual(energy_before, 15.995466129, places=5)
        self.assertAlmostEqual(energy_after, -0.146040876947, places=5)

    def test_ATAT_CG_templates_FAdResS(self):
        # add atomistic interactions
        interNBat = espressopp.interaction.VerletListAdressATLenJonesReacFieldGen(self.vl, self.ftpl)
        potLJ  = espressopp.interaction.LennardJones(epsilon=0.650299305951, sigma=0.316549165245, shift='auto', cutoff=1.4)
        potQQ  = espressopp.interaction.ReactionFieldGeneralized(prefactor=138.935485, kappa=0.0, epsilon1=1.0, epsilon2=80.0, cutoff= 1.4, shift="auto")
        interNBat.setPotential1(type1=0, type2=0, potential=potLJ)
        interNBat.setPotential2(type1=0, type2=0, potential=potQQ)
        self.system.addInteraction(interNBat)

        # add coarse-grained interaction
        interNBcg = espressopp.interaction.VerletListAdressCGHarmonic(self.vl, self.ftpl)
        potCG = espressopp.interaction.Harmonic(K=500.0, r0=1.4, cutoff=1.4)
        interNBcg.setPotential(type1=1, type2=1, potential=potCG)
        self.system.addInteraction(interNBcg)

        # decompose properly
        espressopp.tools.AdressDecomp(self.system, self.integrator)

        # coordinates and non-bonded energy of particles before integration
        before = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_before = interNBcg.computeEnergy() + interNBat.computeEnergy()

        # run ten steps and compute energy
        self.integrator.run(10)

        # coordinates and non-bonded energy of particles after integration
        after = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_after = interNBcg.computeEnergy() + interNBat.computeEnergy()

        # run checks
        self.assertAlmostEqual(after[0], 3.897760, places=5)
        self.assertEqual(before[1], after[1])
        self.assertEqual(before[2], after[2])
        self.assertAlmostEqual(after[3], 4.985460, places=5)
        self.assertEqual(before[4], after[4])
        self.assertEqual(before[5], after[5])
        self.assertAlmostEqual(after[6], 6.947461, places=5)
        self.assertEqual(before[7], after[7])
        self.assertEqual(before[8], after[8])
        self.assertAlmostEqual(after[9], 9.364748, places=5)
        self.assertEqual(before[10], after[10])
        self.assertEqual(before[11], after[11])
        self.assertAlmostEqual(after[12], 2.304571, places=5)
        self.assertEqual(before[13], after[13])
        self.assertEqual(before[14], after[14])
        self.assertAlmostEqual(energy_before, 223.764297, places=5)
        self.assertAlmostEqual(energy_after, 9.190505, places=5)

    def test_AT_CG_templates_FAdResS(self):
        # add atomistic interaction
        interNBat = espressopp.interaction.VerletListAdressATLennardJones(self.vl, self.ftpl)
        potLJ  = espressopp.interaction.LennardJones(epsilon=0.650299305951, sigma=0.316549165245, shift='auto', cutoff=1.4)
        interNBat.setPotential(type1=0, type2=0, potential=potLJ)
        self.system.addInteraction(interNBat)

        # add coarse-grained interaction
        interNBcg = espressopp.interaction.VerletListAdressCGHarmonic(self.vl, self.ftpl)
        potCG = espressopp.interaction.Harmonic(K=500.0, r0=1.4, cutoff=1.4)
        interNBcg.setPotential(type1=1, type2=1, potential=potCG)
        self.system.addInteraction(interNBcg)

        # decompose properly
        espressopp.tools.AdressDecomp(self.system, self.integrator)

        # coordinates and non-bonded energy of particles before integration
        before = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_before = interNBcg.computeEnergy() + interNBat.computeEnergy()

        # run ten steps and compute energy
        self.integrator.run(10)

        # coordinates and non-bonded energy of particles after integration
        after = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_after = interNBcg.computeEnergy() + interNBat.computeEnergy()

        # run checks
        self.assertAlmostEqual(after[0], 9.364194, places=5)
        self.assertEqual(before[1], after[1])
        self.assertEqual(before[2], after[2])
        self.assertAlmostEqual(after[3], 7.260933, places=5)
        self.assertEqual(before[4], after[4])
        self.assertEqual(before[5], after[5])
        self.assertAlmostEqual(after[6], 4.912661, places=5)
        self.assertEqual(before[7], after[7])
        self.assertEqual(before[8], after[8])
        self.assertAlmostEqual(after[9], 0.700571, places=5)
        self.assertEqual(before[10], after[10])
        self.assertEqual(before[11], after[11])
        self.assertAlmostEqual(after[12], 5.261640, places=5)
        self.assertEqual(before[13], after[13])
        self.assertEqual(before[14], after[14])
        self.assertAlmostEqual(energy_before, 199.996600, places=5)
        self.assertAlmostEqual(energy_after, 1.382440, places=5)

if __name__ == '__main__':
    unittest.main()
