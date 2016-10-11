#import sys
#import time
import os
import espressopp
import mpi4py.MPI as MPI
from espressopp import Int3D
from espressopp import Real3D

import unittest

runSteps = 1000
temperature = 1.0
Ni = 10
initDen = 1.
initVel = 0.

class TestPureLB(unittest.TestCase):
    def setUp(self):
        # set up system
        global Ni
        global temperature
        system, integrator = espressopp.standard_system.LennardJones(0, box=(Ni, Ni, Ni), temperature=temperature)
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)

        # set up LB fluid
        lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
        integrator.addExtension(lb)
        lb.lbTemp = 1.0

        # set initial populations
        global initDen, initVel
        initPop = espressopp.integrator.LBInitPopUniform(system,lb)
        initPop.createDenVel(initDen, Real3D(initVel))

        self.system = system
        self.lb = lb
        self.integrator = integrator

    def test_purelb(self):
        global runSteps

        self.integrator.run(runSteps)
        
        # variables to hold average density and mass flux
        av_den = 0.
        av_j = Real3D(0.)
        list = []

        # lattice variables. halo is hard coded
        halo = 1
        myNi = self.lb.getMyNi
        area_yz = (myNi[1] - 2 * halo) * (myNi[2] - 2 * halo)

        for i in range (halo, myNi[0]-halo):
            for j in range (halo, myNi[1]-halo):
                for k in range (halo, myNi[2]-halo):
                    av_den += self.lb.getLBMom(Int3D(i,j,k), 0)

                    jx = self.lb.getLBMom(Int3D(i,j,k), 1)
                    jy = self.lb.getLBMom(Int3D(i,j,k), 2)
                    jz = self.lb.getLBMom(Int3D(i,j,k), 3)
                    av_j += Real3D(jx, jy, jz)
            av_den /= area_yz
            av_j /= area_yz

            self.assertAlmostEqual(av_den, initDen, places=2)
            self.assertAlmostEqual(av_j[0], initVel, places=2)
            self.assertAlmostEqual(av_j[1], initVel, places=2)
            self.assertAlmostEqual(av_j[2], initVel, places=2)

            av_den = 0.
            av_j = Real3D(0.)

if __name__ == '__main__':
    unittest.main()
