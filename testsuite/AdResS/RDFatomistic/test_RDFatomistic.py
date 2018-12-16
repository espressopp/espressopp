#!/usr/bin/env python2
#
#  Copyright (C) 2013-2018
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

class TestRDFatomistic(unittest.TestCase):
    def setUp(self):
        # set up system
        self.system = espressopp.System()
        box = (10, 10, 10)
        self.system.bc = espressopp.bc.OrthorhombicBC(self.system.rng, box)
        self.system.rng = espressopp.esutil.RNG()
        self.system.skin = 0.5
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc=1.5,skin=self.system.skin)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc=1.5, skin=self.system.skin)
        self.system.storage = espressopp.storage.DomainDecompositionAdress(self.system, nodeGrid, cellGrid)

    def test_span_true(self):
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

        # calculate rdfs
        rdf_span = espressopp.analysis.RDFatomistic(system=self.system, type1=0, type2=0, span=1.0, spanbased=True)
        rdf_span_array = rdf_span.compute(10)

        # run checks
        self.assertAlmostEqual(rdf_span_array[1], 0.0, places=5)
        self.assertAlmostEqual(rdf_span_array[2], 33.506304, places=5)

    def test_span_false(self):
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

        # calculate rdfs
        rdf_nospan = espressopp.analysis.RDFatomistic(system=self.system, type1=0, type2=0, span=1.0, spanbased=False)
        rdf_nospan_array = rdf_nospan.compute(10)

        # run checks
        self.assertAlmostEqual(rdf_nospan_array[1], 136.418523, places=5)
        self.assertAlmostEqual(rdf_nospan_array[2], 16.753152, places=5)

    def test_path_integral(self):
        # add some particles
        particle_list = [
            (1, 1, espressopp.Real3D(1.0, 5.0, 5.0), 0, 0),
            (2, 1, espressopp.Real3D(1.0, 4.9, 5.0), 1, 1),
            (3, 1, espressopp.Real3D(1.0, 5.0, 4.9), 1, 2),
            (4, 1, espressopp.Real3D(1.0, 5.1, 5.0), 1, 3),
            (5, 1, espressopp.Real3D(1.0, 5.0, 5.1), 1, 4),
            (6, 1, espressopp.Real3D(6.0, 5.0, 5.0), 0, 0),
            (7, 1, espressopp.Real3D(6.0, 4.8, 5.0), 1, 1),
            (8, 1, espressopp.Real3D(6.0, 5.0, 4.8), 1, 2),
            (9, 1, espressopp.Real3D(6.0, 5.2, 5.0), 1, 3),
            (10, 1, espressopp.Real3D(6.0, 5.0, 5.2), 1, 4),
            (11, 1, espressopp.Real3D(7.0, 5.0, 5.0), 0, 0),
            (12, 1, espressopp.Real3D(7.0, 4.8, 5.0), 1, 1),
            (13, 1, espressopp.Real3D(7.0, 5.0, 4.8), 1, 2),
            (14, 1, espressopp.Real3D(7.0, 5.2, 5.0), 1, 3),
            (15, 1, espressopp.Real3D(7.0, 5.0, 5.2), 1, 4),
            (16, 0, espressopp.Real3D(9.0, 5.0, 5.0), 0, 0),
            (17, 0, espressopp.Real3D(9.0, 4.8, 5.0), 1, 1),
            (18, 0, espressopp.Real3D(9.0, 5.0, 4.8), 1, 2),
            (19, 0, espressopp.Real3D(9.0, 5.2, 5.0), 1, 3),
            (20, 0, espressopp.Real3D(9.0, 5.0, 5.2), 1, 4),
            (21, 0, espressopp.Real3D(5.0, 5.0, 5.0), 0, 0),
            (22, 0, espressopp.Real3D(5.0, 4.8, 5.0), 1, 1),
            (23, 0, espressopp.Real3D(5.0, 5.0, 4.8), 1, 2),
            (24, 0, espressopp.Real3D(5.0, 5.2, 5.0), 1, 3),
            (25, 0, espressopp.Real3D(5.0, 5.0, 5.2), 1, 4),
            (26, 1, espressopp.Real3D(5.0, 5.0, 5.0), 0, 0),
            (27, 1, espressopp.Real3D(5.5, 4.5, 5.0), 1, 2),
            (28, 1, espressopp.Real3D(5.5, 5.0, 4.5), 1, 3),
            (29, 1, espressopp.Real3D(5.5, 5.5, 5.0), 1, 4),
            (30, 1, espressopp.Real3D(5.5, 5.0, 5.5), 1, 1)
        ]
        tuples = [(1,2,3,4,5),(6,7,8,9,10),(11,12,13,14,15),(16,17,18,19,20),(21,22,23,24,25),(26,27,28,29,30)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'pos', 'adrat', 'pib')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)

        # generate a verlet list
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5, dEx=2.0, dHy=1.0, adrCenter=[5.0, 5.0, 5.0], sphereAdr=False)

        # initialize lambda values
        integrator = espressopp.integrator.PIAdressIntegrator(system=self.system, verletlist=vl, nTrotter=4)
        espressopp.tools.AdressDecomp(self.system, integrator)

        # calculate rdfs
        rdf_pi = espressopp.analysis.RDFatomistic(system=self.system, type1=1, type2=1, span=2.5)
        rdf_pi_array = rdf_pi.computePathIntegral(10)

        # run checks
        self.assertAlmostEqual(rdf_pi_array[1], 49.121896, places=5)
        self.assertAlmostEqual(rdf_pi_array[2], 26.525824, places=5)
        self.assertAlmostEqual(rdf_pi_array[3], 15.898631, places=5)
        self.assertAlmostEqual(rdf_pi_array[4], 0.0, places=5)

if __name__ == '__main__':
    unittest.main()
