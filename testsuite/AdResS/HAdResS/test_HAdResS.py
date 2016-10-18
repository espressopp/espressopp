import sys
import time
import espressopp
import mpi4py.MPI as MPI

import unittest

class TestHAdResS(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        box = (10, 10, 10)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, 1.5, 0.3)
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

        # add interaction
        interNB = espressopp.interaction.VerletListHadressLennardJones2(vl, ftpl)
        potWCA1  = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift='auto', cutoff=1.4)
        potWCA2 = espressopp.interaction.LennardJones(epsilon=0.0, sigma=1.0, shift='auto', cutoff=1.4)
        interNB.setPotentialAT(type1=0, type2=0, potential=potWCA1) # AT
        interNB.setPotentialCG(type1=1, type2=1, potential=potWCA2) # CG
        self.system.addInteraction(interNB)

        # initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        # coordinates and non-bonded energy of particles before integration
        before = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_before = interNB.computeEnergy()

        # run ten steps and compute energy
        integrator.run(10)

        # coordinates and non-bonded energy of particles after integration
        after = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_after = interNB.computeEnergy()

        # run checks (only particles in atomistic, hybrid and first part of coarse-grained region should move, i.e. the first 4 particles - as the WCA's epsilon is 0.0 in the coarse-grained region. Furthermore they should move along the x-axis only given their setup. Additionally, check energies)
        self.assertAlmostEqual(after[0], 5.413650, places=5)
        self.assertEqual(before[1], after[1])
        self.assertEqual(before[2], after[2])
        self.assertAlmostEqual(after[3], 6.507075, places=5)
        self.assertEqual(before[4], after[4])
        self.assertEqual(before[5], after[5])
        self.assertAlmostEqual(after[6], 7.551088, places=5)
        self.assertEqual(before[7], after[7])
        self.assertEqual(before[8], after[8])
        self.assertAlmostEqual(after[9], 8.531737, places=5)
        self.assertEqual(before[10], after[10])
        self.assertEqual(before[11], after[11])
        self.assertEqual(before[12], after[12])
        self.assertEqual(before[13], after[13])
        self.assertEqual(before[14], after[14])

        self.assertAlmostEqual(energy_before,0.921374, places=5)
        self.assertAlmostEqual(energy_after, -0.468436, places=5)

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

        # add interaction
        interNB = espressopp.interaction.VerletListHadressLennardJones2(vl, ftpl)
        potWCA1  = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift='auto', cutoff=1.4)
        potWCA2 = espressopp.interaction.LennardJones(epsilon=0.0, sigma=1.0, shift='auto', cutoff=1.4)
        interNB.setPotentialAT(type1=0, type2=0, potential=potWCA1) # AT
        interNB.setPotentialCG(type1=0, type2=0, potential=potWCA2) # CG
        self.system.addInteraction(interNB)

        # initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        # coordinates and non-bonded energy of particles before integration
        before = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_before = interNB.computeEnergy()

        # run ten steps and compute energy
        integrator.run(10)

        # coordinates and non-bonded energy of particles after integration
        after = [self.system.storage.getParticle(i).pos[j] for i in range(1,6) for j in range(3)]
        energy_after = interNB.computeEnergy()

        # run checks (as before, just that particles should now move along the y-axis only, given their setup and the spherical adaptive resolution geometry)
        self.assertEqual(before[0], after[0])
        self.assertAlmostEqual(after[1], 5.413650, places=5)
        self.assertEqual(before[2], after[2])
        self.assertEqual(before[3], after[3])
        self.assertAlmostEqual(after[4], 6.507075, places=5)
        self.assertEqual(before[5], after[5])
        self.assertEqual(before[6], after[6])
        self.assertAlmostEqual(after[7], 7.551088, places=5)
        self.assertEqual(before[8], after[8])
        self.assertEqual(before[9], after[9])
        self.assertAlmostEqual(after[10], 8.531737, places=5)
        self.assertEqual(before[11], after[11])
        self.assertEqual(before[12], after[12])
        self.assertEqual(before[13], after[13])
        self.assertEqual(before[14], after[14])
        self.assertAlmostEqual(energy_before,0.921374, places=5)
        self.assertAlmostEqual(energy_after, -0.468436, places=5)


if __name__ == '__main__':
    unittest.main()
