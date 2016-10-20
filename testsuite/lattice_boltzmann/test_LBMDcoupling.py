import espressopp
import mpi4py.MPI as MPI
from espressopp import Int3D
from espressopp import Real3D

import unittest
import os

runSteps = 500
temperature = 1.0
Npart = 80
Ni = 5
initDen = 1.
initVel = 0.
mdoutput = 'eq_LJ_fluid.xyz'

class makeConf(unittest.TestCase):
    def setUp(self):
        # globals
        global Ni, temperature

        # constants
        rc   = pow(2, 1./6.)
        skin = 0.3
        epsilon = 1.

        # system set up
        system         = espressopp.System()
        system.rng     = espressopp.esutil.RNG()

        global Npart
        if ( os.path.isfile(mdoutput) ):
            pid, ptype, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz(mdoutput)
            Npart = len(pid)
            sigma  = 1.
            timestep = 0.005
            box  = (Lx, Ly, Lz)
        else:
            sigma  = 0.
            timestep = 0.001
            box  = (Ni, Ni, Ni)

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
        integrator.dt  = timestep

        # thermostat
        thermostat     = espressopp.integrator.LangevinThermostat(system)
        thermostat.gamma  = 1.0
        thermostat.temperature = temperature
        integrator.addExtension(thermostat)

        if ( os.path.isfile(mdoutput) ):
            props = ['id', 'type', 'mass', 'pos', 'v']
            new_particles = []
            for pid in range(Npart):
                part = [pid + 1, 0, 1.0, Real3D(x[pid], y[pid], z[pid]), Real3D(vx[pid], vy[pid], vz[pid])]
                new_particles.append(part)
            system.storage.addParticles(new_particles, *props)
            system.storage.decompose()
        else:
            # make dense system
            particle_list = []
            mass = 1.
            for k in range(Npart):
                pid  = k + 1
                pos  = system.bc.getRandomPos()
                v    = Real3D(0.)
                ptype = 0
                part = [pid, ptype, mass, pos, v]
                particle_list.append(part)
            system.storage.addParticles(particle_list, 'id', 'type', 'mass', 'pos', 'v')
            system.storage.decompose()

            print "Warm up. Sigma will be increased from 0. to 1."
            new_sigma = sigma
            for k in range(50):
                integrator.run(100)
                new_sigma += 0.02
                if (new_sigma > 0.8):
                    integrator.dt = 0.005
                potLJ = espressopp.interaction.LennardJones(epsilon, new_sigma, rc)
                interaction.setPotential(type1=0, type2=0, potential=potLJ)

            espressopp.tools.fastwritexyz(mdoutput, system)

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
        ext_lboutputScreen=espressopp.integrator.ExtAnalyze(lboutputScreen,runSteps)
        integrator.addExtension(ext_lboutputScreen)
        
        # lb parameters: viscosities, temperature, time contrast b/w MD and LB
        lb.visc_b = 3.
        lb.visc_s = 3.
        lb.lbTemp = temperature
        lb.nSteps = 5

        # set self
        self.integrator = integrator
        self.lboutput = lboutputScreen

class TestLBMDCoupling(makeConf):
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
