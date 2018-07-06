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


import espressopp
import unittest

Npart = 10
amp = 0.5
vecR = espressopp.Real3D(amp)
 
class Test_CM(unittest.TestCase): 
    def setUp(self):

        system, integrator = espressopp.standard_system.LennardJones(Npart, box=(10,10,10))

        # fill in particles' velocities with dummy numbers
        for i in range (1,Npart + 1):
            system.storage.modifyParticle(i, 'v', i * vecR)
            p = system.storage.getParticle(i)

        # set up analysis of the CMVelocity
        tvel= espressopp.analysis.CMVelocity(system)

        # set self:
        self.system = system
        self.integrator = integrator
        self.tvel = tvel

    def test_initVelocity(self):
        # compute CM-velocity
        self.tvel.compute()

        # check CM-velocity and arithmetic progression's (AP) sum
        ap_sum = Npart * (1 + Npart) * amp / 2.
        self.assertAlmostEqual(self.tvel.v[0], ap_sum / Npart, places=10)
        self.assertAlmostEqual(self.tvel.v[1], ap_sum / Npart, places=10)
        self.assertAlmostEqual(self.tvel.v[2], ap_sum / Npart, places=10)


    def test_CMVelocity(self):
        # attach to integrator
        ext_remove_com = espressopp.integrator.ExtAnalyze(self.tvel, 10)
        self.integrator.addExtension(ext_remove_com)

        for i in range (10):
            self.integrator.run(10)

            # check resetted CM-velocity
            self.assertAlmostEqual(self.tvel.v[0], 0., places=10)
            self.assertAlmostEqual(self.tvel.v[1], 0., places=10)
            self.assertAlmostEqual(self.tvel.v[2], 0., places=10)

            # fill in particles' velocity again with non-zero numbers
            self.system.storage.modifyParticle(i + 1, 'v', (i + 1) * vecR)

if __name__ == '__main__':
    unittest.main()
