#import sys
#import time
import os
import espressopp
import mpi4py.MPI as MPI
from espressopp import Int3D
from espressopp import Real3D

import unittest

runSteps = 500
temperature = 1.0
num_particles = 400
Ni = 10
initDen = 1.
initVel = 0.

class TestLBMDCoupling(unittest.TestCase):
    def setUp(self):
        # set up system
        global Ni, temperature

        box  = (Ni, Ni, Ni)
        rc   = pow(2, 1./6.)
        skin = 0.3
        epsilon = 1.
        sigma  = 0.

        # system set up
        system         = espressopp.System()
        system.rng     = espressopp.esutil.RNG()
        system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin    = skin
        nodeGrid       = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

        # interaction
        interaction    = espressopp.interaction.VerletListLennardJones(espressopp.VerletList(system, cutoff=rc))
        potLJ          = espressopp.interaction.LennardJones(epsilon, sigma, rc)
        interaction.setPotential(type1=0, type2=0, potential=potLJ)
        system.addInteraction(interaction)
        
        # integrator
        integrator     = espressopp.integrator.VelocityVerlet(system)
        integrator.dt = 0.001

        # thermostat
        thermostat     = espressopp.integrator.LangevinThermostat(system)
        thermostat.gamma  = 1.0
        thermostat.temperature = temperature
        integrator.addExtension(thermostat)

        # make dense system
        global num_particles
        particle_list = []
        mass = 1.
        for k in range(num_particles):
            pid  = k + 1
            pos  = system.bc.getRandomPos()
            v    = Real3D(0,0,0)
            type = 0
            part = [pid, pos, type, v, mass]
            particle_list.append(part)
        system.storage.addParticles(particle_list, 'id', 'pos', 'type', 'v', 'mass')
        system.storage.decompose()

        print "Warm up. Sigma will be increased from 0. to 1."
        new_sigma = sigma
        for k in range(100):
            integrator.run(100)
            new_sigma += 0.01
            if (new_sigma > 0.8):
                integrator.dt = 0.005
            potLJ = espressopp.interaction.LennardJones(epsilon, new_sigma, rc)
            interaction.setPotential(type1=0, type2=0, potential=potLJ)

        # LB will control thermostatting now
        thermostat.disconnect()        

        integrator.step = 0

        # set up LB fluid
        lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
        integrator.addExtension(lb)

        # set initial populations
        global initDen, initVel
        initPop = espressopp.integrator.LBInitPopUniform(system,lb)
        initPop.createDenVel(initDen, Real3D(initVel))

        # set up LB profiler
        global runSteps
        lb.profStep = int(.5 * runSteps)
        lboutputScreen = espressopp.analysis.LBOutputScreen(system,lb)
        OUT3=espressopp.integrator.ExtAnalyze(lboutputScreen,runSteps)
        integrator.addExtension(OUT3)
        
        # lb parameters: viscosities, temperature, time contrast b/w MD and LB
        lb.visc_b = 3.
        lb.visc_s = 3.
        lb.lbTemp = temperature
        lb.nSteps = 5

        # set self
        self.lb = lb
        self.lboutput = lboutputScreen
        self.integrator = integrator

    def test_lbmdcoupling(self):
        print "Checking total momentum of LB-MD system:"

        global runSteps

        for checks in range(2):
            self.integrator.run(runSteps)

            # output formatting
            print "-" * 73
            tot_mom = self.lboutput.getLBMom() + self.lboutput.getMDMom()
            print "total LB-MD mom:  ", ("{:>18.1e}"*3).format(*tot_mom), "\n"

            # momentum checks
            self.assertAlmostEqual(tot_mom[0], 0., places=10)
            self.assertAlmostEqual(tot_mom[1], 0., places=10)
            self.assertAlmostEqual(tot_mom[2], 0., places=10)

if __name__ == '__main__':
    unittest.main()
