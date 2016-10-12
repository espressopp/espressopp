import espressopp
import mpi4py.MPI as MPI
from espressopp import Int3D
from espressopp import Real3D

import unittest

runSteps = 1000
temperature = 1.0
Ni = 6
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

        self.lb = lb
        self.integrator = integrator
        self.lbforce = lbforce

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

    def test_extforce(self):
        global runSteps

        vz_check = 0.01
        self.integrator.run(runSteps)
        self.run_average(vz_check)

        # add external constant (gravity-like) force to the existing one
        self.lbforce.addForce(Real3D(0.,0.,0.00002))

        vz_check = 0.04
        self.integrator.run(runSteps)
        self.run_average(vz_check)

if __name__ == '__main__':
    unittest.main()
