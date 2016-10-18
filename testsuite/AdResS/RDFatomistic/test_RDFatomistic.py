import sys
import time
import espressopp
import mpi4py.MPI as MPI

import unittest

class TestRDFatomistic(unittest.TestCase):
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

        # add some particles
        particle_list = [
            (1, 1, 0, espressopp.Real3D(5.0, 5.0, 5.0), 1.0, 0),
            (2, 1, 0, espressopp.Real3D(6.0, 5.0, 5.0), 1.0, 0),
            (3, 1, 0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 0),
            (4, 1, 0, espressopp.Real3D(6.5, 5.5, 5.0), 1.0, 0),
            (5, 0, 0, espressopp.Real3D(5.0, 5.0, 5.0), 1.0, 1),
            (6, 0, 0, espressopp.Real3D(6.0, 5.0, 5.0), 1.0, 1),
            (7, 0, 0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 1),
            (8, 0, 0, espressopp.Real3D(6.5, 5.5, 5.0), 1.0, 1),
        ]
        tuples = [(1,5),(2,6),(3,7),(4,8)]
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
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

    def test_span_true(self):
        # calculate rdfs
        rdf_span = espressopp.analysis.RDFatomistic(system = self.system, type1 = 0, type2 = 0, spanbased = True, span = 1.0)
        rdf_span_array = rdf_span.compute(10)

        # run checks
        self.assertAlmostEqual(rdf_span_array[1], 0.0, places=5)
        self.assertAlmostEqual(rdf_span_array[2], 33.506304, places=5)

    def test_span_false(self):
        # calculate rdfs
        rdf_nospan = espressopp.analysis.RDFatomistic(system = self.system, type1 = 0, type2 = 0, spanbased = False)
        rdf_nospan_array = rdf_nospan.compute(10)

        # run checks
        self.assertAlmostEqual(rdf_nospan_array[1], 136.418523, places=5)
        self.assertAlmostEqual(rdf_nospan_array[2], 16.753152, places=5)


if __name__ == '__main__':
    unittest.main()
