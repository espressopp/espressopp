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


class TestThermoIntegration(unittest.TestCase):
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

    def test_potential_decoupling(self):
        #add particles
        particle_list = [
            (4, 1,  0, espressopp.Real3D(2.0, 2.0, 2.0), 1.0, 0),
            (5, 1,  0, espressopp.Real3D(2.3, 2.0, 2.0), 1.0, 0),
            (6, 1,  0, espressopp.Real3D(2.6, 2.0, 2.0), 1.0, 0),
            (1, 0,  1, espressopp.Real3D(2.0, 2.0, 2.0), 1.0, 1),
            (2, 0, -1, espressopp.Real3D(2.3, 2.0, 2.0), 1.0, 1),
            (3, 0, -1, espressopp.Real3D(2.6, 2.0, 2.0), 1.0, 1),
        ]
        tuples = [(4,1),(5,2),(6,3)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                dEx=2.0, dHy=1.0, pids=[4], sphereAdr=True)

        #add LJ interaction
        interactionLJ = espressopp.interaction.VerletListAdressLennardJonesSoftcoreTI(vl, ftpl)
        potLJ = espressopp.interaction.LennardJonesSoftcoreTI(epsilonA=1.0, sigmaA=0.2, epsilonB=0.0, sigmaB=0.2, alpha=0.5, power=1.0, cutoff=1.5, lambdaTI=0.3, annihilate=False)
        potLJ.addPids([1,2])
        interactionLJ.setPotentialAT(type1=0, type2=0, potential=potLJ)
        self.system.addInteraction(interactionLJ)

        #add electrostatic interaction
        interactionQQ = espressopp.interaction.VerletListAdressReactionFieldGeneralizedTI(vl, ftpl)
        potQQ = espressopp.interaction.ReactionFieldGeneralizedTI(prefactor=1.0, kappa=0.0, epsilon1=1, epsilon2=80, cutoff=1.5, lambdaTI=0.3, annihilate=False)
        potQQ.addPids([1,2])
        interactionQQ.setPotentialAT(type1=0, type2=0, potential=potQQ)
        self.system.addInteraction(interactionQQ)

        #initialize lambda values
        integrator     = espressopp.integrator.VelocityVerlet(self.system)
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        self.assertAlmostEqual(interactionQQ.computeEnergyDeriv(),-1.627412,places=5)
        self.assertAlmostEqual(interactionLJ.computeEnergyDeriv(),0.330739,places=5)
        self.assertAlmostEqual(interactionQQ.computeEnergy(),-1.213441,places=5)
        self.assertAlmostEqual(interactionLJ.computeEnergy(),-0.545769,places=5)

    def test_potential_annihilation(self):
        #add particles
        particle_list = [
            (4, 1,  0, espressopp.Real3D(2.0, 2.0, 2.0), 1.0, 0),
            (5, 1,  0, espressopp.Real3D(2.3, 2.0, 2.0), 1.0, 0),
            (6, 1,  0, espressopp.Real3D(2.6, 2.0, 2.0), 1.0, 0),
            (1, 0,  1, espressopp.Real3D(2.0, 2.0, 2.0), 1.0, 1),
            (2, 0, -1, espressopp.Real3D(2.3, 2.0, 2.0), 1.0, 1),
            (3, 0, -1, espressopp.Real3D(2.6, 2.0, 2.0), 1.0, 1),
        ]
        tuples = [(4,1),(5,2),(6,3)]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat')
        ftpl = espressopp.FixedTupleListAdress(self.system.storage)
        ftpl.addTuples(tuples)
        self.system.storage.setFixedTuplesAdress(ftpl)
        self.system.storage.decompose()
        vl = espressopp.VerletListAdress(self.system, cutoff=1.5, adrcut=1.5,
                                dEx=2.0, dHy=1.0, pids=[4], sphereAdr=True)

        #add LJ interaction
        interactionLJ = espressopp.interaction.VerletListAdressLennardJonesSoftcoreTI(vl, ftpl)
        potLJ = espressopp.interaction.LennardJonesSoftcoreTI(epsilonA=1.0, sigmaA=0.2, epsilonB=0.0, sigmaB=0.2, alpha=0.5, power=1.0, cutoff=1.5, lambdaTI=0.3, annihilate=True)
        potLJ.addPids([1,2])
        interactionLJ.setPotentialAT(type1=0, type2=0, potential=potLJ)
        self.system.addInteraction(interactionLJ)

        #add electrostatic interaction
        interactionQQ = espressopp.interaction.VerletListAdressReactionFieldGeneralizedTI(vl, ftpl)
        potQQ = espressopp.interaction.ReactionFieldGeneralizedTI(prefactor=1.0, kappa=0.0, epsilon1=1, epsilon2=80, cutoff=1.5, lambdaTI=0.3, annihilate=True)
        potQQ.addPids([1,2])
        interactionQQ.setPotentialAT(type1=0, type2=0, potential=potQQ)
        self.system.addInteraction(interactionQQ)

        #initialize lambda values
        integrator     = espressopp.integrator.VelocityVerlet(self.system)
        adress = espressopp.integrator.Adress(self.system,vl,ftpl)
        integrator.addExtension(adress)
        espressopp.tools.AdressDecomp(self.system, integrator)

        self.assertAlmostEqual(interactionQQ.computeEnergyDeriv(),0.725217,places=5)
        self.assertAlmostEqual(interactionLJ.computeEnergyDeriv(),0.655998,places=5)
        self.assertAlmostEqual(interactionQQ.computeEnergy(),-0.507652,places=5)
        self.assertAlmostEqual(interactionLJ.computeEnergy(),-0.447031,places=5)


if __name__ == '__main__':
    unittest.main()
