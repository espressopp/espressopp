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

class TestCaseConstrainCOM(unittest.TestCase):

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
	self.skin =system.skin

    def test_constrain_com(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,self.box,rc=1.5, skin=self.skin)
        cellGrid = espressopp.tools.decomp.cellGrid(self.box, nodeGrid, rc=1.5, skin=self.skin)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        particle_list = [
            (1, 1, 0, espressopp.Real3D(3.0, 5.0, 5.0), 1.0, 0, 1.),
            (2, 1, 0, espressopp.Real3D(4.0, 5.0, 5.0), 2.0, 0, 1.),
            (3, 1, 0, espressopp.Real3D(5.0, 5.0, 5.0), 3.0, 0, 1.),
            (4, 1, 0, espressopp.Real3D(6.0, 5.0, 5.0), 2.0, 0, 1.),
            (5, 1, 0, espressopp.Real3D(7.0, 5.0, 5.0), 1.0, 0, 1.),
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
        potCOM = espressopp.interaction.ConstrainCOM(1000.)
        interCOM = espressopp.interaction.FixedLocalTupleListConstrainCOM(self.system, tuplelist, potCOM)
        self.system.addInteraction(interCOM, 'Constrain_COM')

        # center of mass of particles before integration
        before = [0., 0., 0.]
        
        particle = self.system.storage.getParticle(1)
        dmy_p = []
        dmy_ele = []
        mass = []
        for i in xrange(3):
            dmy_ele.append(particle.pos[i])
        dmy_p.append(dmy_ele)
        mass.append(particle.mass)
        for i in xrange(2, 6):
            particle = self.system.storage.getParticle(i)
            mass.append(particle.mass)
            diff = []
            for j in xrange(3):
                x_i = particle.pos[j] - dmy_p[i - 2][j]
                x_i = x_i - round(x_i/self.L)*self.L
                diff.append(x_i + dmy_p[i - 2][j])
            dmy_p.append(diff)
        total_mass = 0.
        for i in xrange(5):
            total_mass += mass[i]
            for j in xrange(3):
                before[j] += mass[i]*dmy_p[i][j]
        for i in xrange(3):
            before[i] /= total_mass
        print "before", before

        # run twenty thousand steps
        integrator.run(20000)

        # center of mass of particles after integration
        after= [0., 0., 0.]

        particle = self.system.storage.getParticle(1)
        dmy_p = []
        dmy_ele = []
        mass = []
        for i in xrange(3):
            dmy_ele.append(particle.pos[i])
        dmy_p.append(dmy_ele)
        mass.append(particle.mass)
        for i in xrange(2, 6):
            particle = self.system.storage.getParticle(i)
            mass.append(particle.mass)
            diff = []
            for j in xrange(3):
                x_i = particle.pos[j] - dmy_p[i - 2][j]
                x_i = x_i - round(x_i/self.L)*self.L
                diff.append(x_i + dmy_p[i - 2][j])
            dmy_p.append(diff)
        total_mass = 0.
        for i in xrange(5):
            total_mass += mass[i]
            for j in xrange(3):
                after[j] += mass[i]*dmy_p[i][j]
        for i in xrange(3):
            after[i] /= total_mass
        print "after", after

        # run checks
        self.assertTrue(fabs(before[0] - after[0]) < 0.04)
        self.assertTrue(fabs(before[1] - after[1]) < 0.04)
        self.assertTrue(fabs(before[2] - after[2]) < 0.04)


if __name__ == '__main__':
    unittest.main()
