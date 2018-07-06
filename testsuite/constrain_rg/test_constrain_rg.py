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

from math import fabs
import sys
import time
import espressopp
import mpi4py.MPI as MPI

import unittest

class TestCaseConstrainRG(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        rng = espressopp.esutil.RNG()
        rng.seed(1)
        system.rng = rng
        L = 10
        box = (L, L, L)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        self.system = system
        self.L = L
        self.box = box

    def test_constrain_rg(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,self.box, rc=1.5, skin=0.3)
        cellGrid = espressopp.tools.decomp.cellGrid(self.box, nodeGrid, rc=1.5, skin=0.3)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        particle_list = [
            (1, 1, 0, espressopp.Real3D(4.0, 5.0, 5.0), 1.0, 0, 1.),
            (2, 1, 0, espressopp.Real3D(4.5, 5.866, 5.0), 1.0, 0, 1.),
            (3, 1, 0, espressopp.Real3D(5.0, 5.0, 5.0), 1.0, 0, 1.),
            (4, 1, 0, espressopp.Real3D(5.5, 4.134, 5.0), 1.0, 0, 1.),
            (5, 1, 0, espressopp.Real3D(6.0, 5.0, 5.0), 1.0, 0, 1.),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat', 'radius')
        self.system.storage.decompose()

        # integrator
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.005

        # Langevin Thermostat
        langevin = espressopp.integrator.LangevinThermostat(self.system)
        langevin.gamma = 1.0
        langevin.temperature = 1.0
        integrator.addExtension(langevin)

        # Harmonic bonds
        bondlist  = espressopp.FixedPairList(self.system.storage)
        for i in range(1, 5):
            bondlist.add(i, i + 1)
        potBond = espressopp.interaction.Harmonic(K=100., r0 = 1.0)
        interBond = espressopp.interaction.FixedPairListHarmonic(self.system, bondlist, potBond)
        self.system.addInteraction(interBond)

        # constrain center of mass
        tuplelist = espressopp.FixedLocalTupleList(self.system.storage)
        tuple = []
        for i in range(1, 6):
            tuple.append(i)
        tuplelist.addTuple(tuple)
        potRG = espressopp.interaction.ConstrainRG(2000.)
        interRG = espressopp.interaction.FixedLocalTupleListConstrainRG(self.system, tuplelist, potRG)
        self.system.addInteraction(interRG, 'Constrain_RG')

        # radius of gyration of particles before integration
        before = [0., 0., 0.]
        
        particle = self.system.storage.getParticle(1)
        dmy_p = []
        dmy_ele = []
        for i in xrange(3):
            dmy_ele.append(particle.pos[i])
        dmy_p.append(dmy_ele)
        for i in xrange(2, 6):
            particle = self.system.storage.getParticle(i)
            diff = []
            for j in xrange(3):
                x_i = particle.pos[j] - dmy_p[i - 2][j]
                x_i = x_i - round(x_i/self.L)*self.L
                diff.append(x_i + dmy_p[i - 2][j])
            dmy_p.append(diff)
        for i in xrange(5):
            for j in xrange(3):
                before[j] += dmy_p[i][j]
        for i in xrange(3):
            before[i] /= 5.
        print "before COM =", before
        before_rg = 0.
        for i in xrange(5):
            for j in xrange(3):
                before_rg += (dmy_p[i][j] - before[j])**2
        before_rg = before_rg**0.5
        print "before Rg =", before_rg

        # run twenty thousand steps
        integrator.run(20000)

        # center of mass of particles after integration
        after= [0., 0., 0.]

        particle = self.system.storage.getParticle(1)
        dmy_p = []
        dmy_ele = []
        for i in xrange(3):
            dmy_ele.append(particle.pos[i])
        dmy_p.append(dmy_ele)
        for i in xrange(2, 6):
            particle = self.system.storage.getParticle(i)
            diff = []
            for j in xrange(3):
                x_i = particle.pos[j] - dmy_p[i - 2][j]
                x_i = x_i - round(x_i/self.L)*self.L
                diff.append(x_i + dmy_p[i - 2][j])
            dmy_p.append(diff)
        for i in xrange(5):
            for j in xrange(3):
                after[j] += dmy_p[i][j]
        for i in xrange(3):
            after[i] /= 5.
        print "after  COM =", after
        after_rg = 0.
        for i in xrange(5):
            for j in xrange(3):
                after_rg += (dmy_p[i][j] - after[j])**2
        after_rg = after_rg**0.5
        print "after  Rg =", after_rg

        # run checks
        self.assertTrue(fabs((before_rg - after_rg)/before_rg) < 0.03)

if __name__ == '__main__':
    unittest.main()
