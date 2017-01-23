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
initVelSin = 0.1

class TestPureLB(unittest.TestCase):
    def setUp(self):
        # set up system
        global Ni, temperature
        system, integrator = espressopp.standard_system.LennardJones(0, box=(Ni, Ni, Ni), temperature=temperature)
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)

        # set up LB fluid
        lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
        integrator.addExtension(lb)

        self.system = system
        self.lb = lb
        self.integrator = integrator

    def test_purelb(self):
        print "Checking athermal LB fluid"
        global runSteps

        # set initial populations
        global initDen, initVel
        initPop = espressopp.integrator.LBInitPopUniform(self.system,self.lb)
        initPop.createDenVel(initDen, Real3D(initVel))

        self.integrator.run(runSteps)

        self.check_averages(initVel)

    def test_thermallb(self):
        print "Checking stochastic LB fluid with sin-wave initialization"
        global runSteps

        # set initial populations
        global initDen, initVel, initVelSin
        initPop = espressopp.integrator.LBInitPopWave(self.system,self.lb)
        initPop.createDenVel(initDen, Real3D(initVel, initVel, initVelSin))

        # set temperature 
        self.lb.lbTemp = temperature
        self.integrator.run(runSteps)

        self.check_averages(initVel) # sin-like wave is killed by temperature

    def check_averages(self, _v):
        # variables to hold average density and mass flux
        av_den = 0.
        av_j = Real3D(0.)

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

            print av_den, av_j

            self.assertAlmostEqual(av_den, initDen, places=2)
            self.assertAlmostEqual(av_j[0], _v, places=2)
            self.assertAlmostEqual(av_j[1], _v, places=2)
            self.assertAlmostEqual(av_j[2], _v, places=2)

            av_den = 0.
            av_j = Real3D(0.)

if __name__ == '__main__':
    unittest.main()
