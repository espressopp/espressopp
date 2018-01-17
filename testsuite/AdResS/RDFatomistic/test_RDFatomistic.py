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

class TestRDFatomistic(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        box = (10, 10, 10)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc=1.5,system.skin)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc=1.5, system.skin)
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
