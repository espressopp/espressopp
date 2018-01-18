#!/usr/bin/env python
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
import math
import unittest


def sqrlen(vector):
    return vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]


class TestRattle(unittest.TestCase):
    def setUp(self):

        system = espressopp.System()
        box = (10, 10, 10)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc=1.5,skin=system.skin)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc=1.5, skin=system.skin)
        system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)
        self.system = system


    def test_rattle(self):

        #add particles
        # 4-3-5
        #   |
        #   |
        # 2-1
        particle_list = [
            (6, 1, espressopp.Real3D(3.2, 3.3, 3.0), espressopp.Real3D(0,0,0),  9.0, 0),
            (1, 0, espressopp.Real3D(3.2, 3.1, 3.0), espressopp.Real3D(0.27, -0.07,  0.10),  3.0, 1),
            (2, 0, espressopp.Real3D(3.1, 3.1, 3.0), espressopp.Real3D(0.23,  0.22,  0.00),  1.0, 1),
            (3, 0, espressopp.Real3D(3.2, 3.3, 3.0), espressopp.Real3D(0.24, -0.37, -0.10),  3.0, 1),
            (4, 0, espressopp.Real3D(3.1, 3.3, 3.0), espressopp.Real3D(0.14,  0.20, -0.05),  1.0, 1),
            (5, 0, espressopp.Real3D(3.3, 3.3, 3.0), espressopp.Real3D(0.04,  0.17, -0.20),  1.0, 1) 
        ]

        tuples = [(6,1,2,3,4,5)]

        constrainedBondsList = [[1, 2, 0.1, 3.0, 1.0], [3, 4, 0.1, 3.0, 1.0], [3, 5, 0.1, 3.0, 1.0]] #pid1,pid2,constraintDist,mass1,mass2
        constraintDist2 = []
        for bond in constrainedBondsList:
          constraintDist2.append(bond[2]*bond[2])

        self.system.storage.addParticles(particle_list, 'id', 'type', 'pos', 'v', 'mass','adrat')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                dEx=2.0, dHy=1.0, pids=[6], sphereAdr=True)

        #one unconstrained bond
        fpl = espressopp.FixedPairListAdress(self.system.storage, ftpl)
        fpl.addBonds([(1,3)])
        pot = espressopp.interaction.Harmonic(K=5.0, r0=0.2) 
        interB = espressopp.interaction.FixedPairListHarmonic(self.system, fpl, pot)
        self.system.addInteraction(interB)
        #some angles
        ftl = espressopp.FixedTripleListAdress(self.system.storage, ftpl)
        ftl.addTriples([(2,1,3),(1,3,4),(1,3,5)])
        pot = espressopp.interaction.AngularHarmonic(K=5.0, theta0=math.pi/2.0) 
        interA = espressopp.interaction.FixedTripleListAngularHarmonic(self.system, ftl, pot)
        self.system.addInteraction(interA)

        #initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        rattle = espressopp.integrator.Rattle(self.system, maxit = 1000, tol = 1e-6, rptol = 1e-6)
        rattle.addConstrainedBonds(constrainedBondsList)
        integrator.addExtension(rattle)

        integrator.run(5)

        #check if bond lengths are the same as constraint lengths
        for i,bond in enumerate(constrainedBondsList):
          pid1 = bond[0]
          pid2 = bond[1]
          dist2 = sqrlen(self.system.bc.getMinimumImageVector(self.system.storage.getParticle(pid1).pos,self.system.storage.getParticle(pid2).pos))
          self.assertAlmostEqual(constraintDist2[i],dist2,places=6)

        #check velocities after 5 steps of deterministic simulation
        vsum = 0.0
        for pid in xrange(1,6):
          part = self.system.storage.getParticle(pid)
          vsum += sqrlen(part.v)
        self.assertAlmostEqual(vsum,0.3842668659,places=6)


if __name__ == '__main__':
    unittest.main()
