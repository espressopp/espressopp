import sys
import time
import espressopp
import mpi4py.MPI as MPI

import unittest

class TestThermoIntegration(unittest.TestCase):
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
                                dEx=1.0, dHy=1.0, adrCenter=[5.0, 5.0, 5.0], sphereAdr=False)

        # initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        # set up TD force
        thdforce = espressopp.integrator.TDforce(self.system,vl)
        thdforce.addForce(itype=3,filename="table_tf.tab",type=1)
        integrator.addExtension(thdforce)

        # x coordinates of particles before integration
        before = [self.system.storage.getParticle(i).pos[0] for i in range(1,6)]

        # run ten steps
        integrator.run(10)

        # x coordinates of particles after integration
        after = [self.system.storage.getParticle(i).pos[0] for i in range(1,6)]

        # run checks
        self.assertAlmostEqual(before[0], after[0], places=5)
        self.assertAlmostEqual(after[1], 6.596913, places=5)
        self.assertAlmostEqual(before[2], after[2], places=5)
        self.assertAlmostEqual(before[3], after[3], places=5)
        self.assertAlmostEqual(before[4], after[4], places=5)

    def test_fixed_sphere(self):
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
                                dEx=1.0, dHy=1.0, adrCenter=[5.0, 5.0, 5.0], sphereAdr=True)

        # initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        # set up TD force
        thdforce = espressopp.integrator.TDforce(self.system,vl)
        thdforce.addForce(itype=3,filename="table_tf.tab",type=1)
        integrator.addExtension(thdforce)

        # y coordinates of particles before integration
        before = [self.system.storage.getParticle(i).pos[1] for i in range(1,6)]

        # run ten steps
        integrator.run(10)

        # y coordinates of particles after integration
        after = [self.system.storage.getParticle(i).pos[1] for i in range(1,6)]

        # run checks
        self.assertAlmostEqual(before[0], after[0], places=5)
        self.assertAlmostEqual(after[1], 6.596913, places=5)
        self.assertAlmostEqual(before[2], after[2], places=5)
        self.assertAlmostEqual(before[3], after[3], places=5)
        self.assertAlmostEqual(before[4], after[4], places=5)

    def test_single_moving_sphere(self):
        # add some particles
        particle_list = [
            (1, 1, 0, espressopp.Real3D(5.0, 5.0, 5.0), espressopp.Real3D(0.0, 1.0, 0.0), 1.0, 0),
            (2, 1, 0, espressopp.Real3D(5.0, 6.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 0),
            (3, 1, 0, espressopp.Real3D(5.0, 7.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 0),
            (4, 1, 0, espressopp.Real3D(5.0, 8.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 0),
            (5, 1, 0, espressopp.Real3D(5.0, 9.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 0),
            (6, 0, 0, espressopp.Real3D(5.0, 5.0, 5.0), espressopp.Real3D(0.0, 1.0, 0.0), 1.0, 1),
            (7, 0, 0, espressopp.Real3D(5.0, 6.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 1),
            (8, 0, 0, espressopp.Real3D(5.0, 7.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 1),
            (9, 0, 0, espressopp.Real3D(5.0, 8.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 1),
            (10, 0, 0, espressopp.Real3D(5.0, 9.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 1),
        ]
        tuples = [(1,6),(2,7),(3,8),(4,9),(5,10)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'v', 'mass','adrat')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()

        # generate a verlet list
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                dEx=1.0, dHy=1.0, pids=[1], sphereAdr=True)

        # initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        # set up TD force
        thdforce = espressopp.integrator.TDforce(self.system,vl)
        thdforce.addForce(itype=3,filename="table_tf.tab",type=1)
        integrator.addExtension(thdforce)

        # y coordinates of particles before integration
        before = [self.system.storage.getParticle(i).pos[1] for i in range(1,6)]

        # run ten steps
        integrator.run(10)

        # y coordinates of particles after integration
        after = [self.system.storage.getParticle(i).pos[1] for i in range(1,6)]

        # run checks
        self.assertAlmostEqual(after[0], 5.100000, places=5)
        self.assertAlmostEqual(after[1], 6.618377, places=5)
        self.assertAlmostEqual(before[2], after[2], places=5)
        self.assertAlmostEqual(before[3], after[3], places=5)
        self.assertAlmostEqual(before[4], after[4], places=5)

    def test_many_moving_spheres(self):
        # add some particles
        particle_list = [
            (1, 1, 0, espressopp.Real3D(4.0, 5.0, 5.0), espressopp.Real3D(-1.0, 1.0, 0.0), 1.0, 0),
            (2, 1, 0, espressopp.Real3D(6.0, 5.0, 5.0), espressopp.Real3D(1.0, 1.0, 0.0), 1.0, 0),
            (3, 1, 0, espressopp.Real3D(5.0, 6.2, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 0),
            (4, 1, 0, espressopp.Real3D(3.0, 8.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 0),
            (5, 1, 0, espressopp.Real3D(7.0, 9.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 0),
            (6, 0, 0, espressopp.Real3D(4.0, 5.0, 5.0), espressopp.Real3D(-1.0, 1.0, 0.0), 1.0, 1),
            (7, 0, 0, espressopp.Real3D(6.0, 5.0, 5.0), espressopp.Real3D(1.0, 1.0, 0.0), 1.0, 1),
            (8, 0, 0, espressopp.Real3D(5.0, 6.2, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 1),
            (9, 0, 0, espressopp.Real3D(3.0, 8.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 1),
            (10, 0, 0, espressopp.Real3D(7.0, 9.5, 5.0), espressopp.Real3D(0.0, 0.0, 0.0), 1.0, 1),
        ]
        tuples = [(1,6),(2,7),(3,8),(4,9),(5,10)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'v', 'mass','adrat')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()

        # generate a verlet list
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                dEx=1.0, dHy=1.0, pids=[1,2], sphereAdr=True)

        # initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        # set up TD force
        thdforce = espressopp.integrator.TDforce(self.system, vl, startdist = 0.9, enddist = 2.1, edgeweightmultiplier = 20)
        thdforce.addForce(itype=3,filename="table_tf.tab",type=1)
        integrator.addExtension(thdforce)

        # y coordinates of particles before integration
        before = [self.system.storage.getParticle(i).pos[1] for i in range(1,6)]

        # run ten steps
        integrator.run(10)

        # y coordinates of particles after integration
        after = [self.system.storage.getParticle(i).pos[1] for i in range(1,6)]

        # run checks
        self.assertAlmostEqual(after[0], 5.100000, places=5)
        self.assertAlmostEqual(after[1], 5.100000, places=5)
        self.assertAlmostEqual(after[2], 6.260323, places=5)
        self.assertAlmostEqual(before[3], after[3], places=5)
        self.assertAlmostEqual(before[4], after[4], places=5)


if __name__ == '__main__':
    unittest.main()
