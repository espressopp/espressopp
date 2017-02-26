#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import time
import espressopp
import mpi4py.MPI as MPI
import math
import unittest


class TestFixed(unittest.TestCase):
    def setUp(self):

        system = espressopp.System()
        box = (10.0, 10.0, 10.0)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        self.cutoff = 4.5
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, self.cutoff, system.skin)
        system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)
        self.system = system


    def test_isfixed(self):

        #add particles
        # AT particles 1 2, VP particle 6 (always atomistic)
        # AT particles 3 4, VP particle 7 (adaptive resolution, based on distance from 5/8)
        # AT particle 5, VP particle 8 (always coarse-grained)

        properties   = ['id', 'type', 'pos', 'v', 'mass', 'adrat', 'isfixed']
        # adrat: 0 is VP-level particle, 1 is AT-level particle
        # isFixed: -1 is adres, 0 is fixed cg, 1 is fixed aa, only value of isFixed for VP-level particle is used, value for corresponding AT-level particle(s) not used

        particle_list = [
            (6, 3, espressopp.Real3D(5.0, 5.0, 5.5), espressopp.Real3D(0,0,0),  2.0, 0,  1),
            (1, 0, espressopp.Real3D(5.0, 5.0, 5.0), espressopp.Real3D(0,0,0),  1.0, 1,  1),
            (2, 1, espressopp.Real3D(5.0, 5.0, 6.0), espressopp.Real3D(0,0,0),  1.0, 1,  1),
            (7, 3, espressopp.Real3D(4.0, 4.0, 3.5), espressopp.Real3D(0,0,0),  2.0, 0, -1),
            (3, 2, espressopp.Real3D(4.0, 4.0, 3.0), espressopp.Real3D(0,0,0),  1.0, 1, -1),
            (4, 0, espressopp.Real3D(4.0, 4.0, 4.0), espressopp.Real3D(0,0,0),  1.0, 1, -1),
            (8, 3, espressopp.Real3D(3.0, 3.0, 3.0), espressopp.Real3D(0,0,0),  1.0, 0,  0),
            (5, 1, espressopp.Real3D(3.0, 3.0, 3.0), espressopp.Real3D(0,0,0),  1.0, 1,  0) 
        ]
        #particle 8 is at centre of AT region
        #distance 8-7 = 1.5
        #distance 8-6 = 3.77

        tuples = [(6,1,2),(7,3,4),(8,5)]

        self.system.storage.addParticles(particle_list, *properties)
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()
        vl = espressopp.VerletListAdress(self.system, cutoff=self.cutoff, adrcut=self.cutoff,
                                dEx=1.0, dHy=1.0, pids=[8], sphereAdr=True)

        #add LJ interaction
        interaction = espressopp.interaction.VerletListAdressLennardJones2(vl, ftpl)
        sigma = [2.4,2.5,2.7]
        epsilon = [1.0,2.0,3.0]
        sigCG = 2.0
        epsCG = 1.0
        for i in xrange(3):
            for j in xrange(3):
                sig = 0.5*(sigma[i]+sigma[j])
                eps = math.sqrt(epsilon[i]*epsilon[j])
                ljpot = espressopp.interaction.LennardJones(epsilon=eps, sigma=sig, shift='auto', cutoff=self.cutoff)
		interaction.setPotentialAT(type1=i, type2=j, potential=ljpot)
        ljpot = espressopp.interaction.LennardJones(epsilon=epsCG, sigma=sigCG, shift='auto', cutoff=self.cutoff)
        interaction.setPotentialCG(type1=3, type2=3, potential=ljpot)

        self.system.addInteraction(interaction)

        #initialize lambda values
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)
        integrator.run(0)

        integrator.run(5)

        self.assertAlmostEqual(interaction.computeEnergy(),21.1042220069,places=6)

        #fixed resolution AT particle in CG region
        part = self.system.storage.getParticle(6)
        self.assertAlmostEqual(part.lambda_adr,1.000000,places=6)

        #adaptive resolution particle in HY region
        part = self.system.storage.getParticle(7)
        self.assertAlmostEqual(part.lambda_adr,0.15923008008,places=6)

        #fixed resolution CG particle in AT region
        part = self.system.storage.getParticle(8)
        self.assertAlmostEqual(part.lambda_adr,0.000000,places=6)


if __name__ == '__main__':
    unittest.main()
