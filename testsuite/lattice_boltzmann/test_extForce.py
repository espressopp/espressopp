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
import mpi4py.MPI as MPI
from espressopp import Int3D
from espressopp import Real3D

import unittest

runSteps = 500
temperature = 1.0
Ni = 5
initDen = 1.
initVel = 0.

class TestExtForceLB(unittest.TestCase):
    def setUp(self):
        # set up system
        global Ni, temperature
        system, integrator = espressopp.standard_system.LennardJones(0, box=(Ni, Ni, Ni), temperature=temperature)
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)

        # set up LB fluid
        lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
        integrator.addExtension(lb)
        lb.lbTemp = temperature

        # set initial populations
        global initDen, initVel
        initPop = espressopp.integrator.LBInitPopUniform(system,lb)
        initPop.createDenVel(initDen, Real3D(initVel))

        # set external constant (gravity-like) force
        lbforce = espressopp.integrator.LBInitConstForce(system,lb)
        lbforce.setForce(Real3D(0.,0.,0.00001))

        # set external sin-like force
        lbforceSin = espressopp.integrator.LBInitPeriodicForce(system,lb)

        self.lb = lb
        self.integrator = integrator
        self.lbforce = lbforce
        self.lbforceSin = lbforceSin

    def run_average(self, _vz):
        # variables to hold average density and mass flux
        av_den = 0.
        av_j = Real3D(0.)

        # lattice variables. halo is hard coded
        halo = 1
        myNi = self.lb.getMyNi
        volume_xyz = (myNi[0] - 2 * halo) * (myNi[1] - 2 * halo) * (myNi[2] - 2 * halo)
        
        for i in range (halo, myNi[0]-halo):
            for j in range (halo, myNi[1]-halo):
                for k in range (halo, myNi[2]-halo):
                    av_den += self.lb.getLBMom(Int3D(i,j,k), 0)

                    jx = self.lb.getLBMom(Int3D(i,j,k), 1)
                    jy = self.lb.getLBMom(Int3D(i,j,k), 2)
                    jz = self.lb.getLBMom(Int3D(i,j,k), 3)
                    av_j += Real3D(jx, jy, jz)
        av_den /= volume_xyz
        av_j /= volume_xyz

        print av_den, av_j

        self.assertAlmostEqual(av_den, initDen, places=2)
        self.assertAlmostEqual(av_j[0], initVel, places=2)
        self.assertAlmostEqual(av_j[1], initVel, places=2)
        self.assertAlmostEqual(av_j[2], _vz, places=2)

    def test_constextforce(self):
        global runSteps

        vz_check = 0.005
        self.integrator.run(runSteps)
        self.run_average(vz_check)

        # add external constant (gravity-like) force to the existing one
        self.lbforce.addForce(Real3D(0.,0.,0.00002))

        vz_check = 0.02
        self.integrator.run(runSteps)
        self.run_average(vz_check)

    def test_sinextforce(self):
        global runSteps
        self.lbforce.addForce(Real3D(0.))

        self.lbforceSin.setForce(Real3D(0.,0.,0.00001))

        vz_check = 0.   # average vz always gives 0 with sin-line force fz(x)
        self.integrator.run(runSteps)
        self.run_average(vz_check)

        # add external constant (gravity-like) force to the existing one
        self.lbforceSin.addForce(Real3D(0.,0.,0.00002))

        self.integrator.run(runSteps)
        self.run_average(vz_check)

if __name__ == '__main__':
    unittest.main()
