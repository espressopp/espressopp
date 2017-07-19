import sys
import time
import espressopp
import mpi4py.MPI as MPI

import unittest

class TestCaseLangevinThermostatOnRadius(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        rng = espressopp.esutil.RNG()
        rng.seed(1)
        system.rng = rng
        box = (10, 10, 10)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        self.system = system

    def test_on_radius(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid((10, 10, 10), nodeGrid, 1.5, 0.3)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        particle_list = [
            (1, 1, 0, espressopp.Real3D(5.5, 5.0, 5.0), 1.0, 0, 1.),
            (2, 1, 0, espressopp.Real3D(6.5, 5.0, 5.0), 1.0, 0, 1.),
            (3, 1, 0, espressopp.Real3D(7.5, 5.0, 5.0), 1.0, 0, 1.),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass','adrat', 'radius')
        self.system.storage.decompose()

        # generate a verlet list
        vl = espressopp.VerletList(self.system, cutoff=1.5)

        # integrator
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01

        # integrator on radius
        radius_mass = 10.
        integratorOnRadius  = espressopp.integrator.VelocityVerletOnRadius(self.system, dampingmass=radius_mass)
        integrator.addExtension(integratorOnRadius)

        # Langevin Thermostat on Radius
        langevin = espressopp.integrator.LangevinThermostatOnRadius(self.system, dampingmass=radius_mass)
        langevin.gamma = 1.0
        langevin.temperature = 10.0
        langevin.adress = False
        langevin.addExclusions([1])
        integrator.addExtension(langevin)

        # coordinates of particles before integration
        before =[self.system.storage.getParticle(i).radius for i in range(1,4)]

        # run ten steps
        integrator.run(10)

        # coordinates of particles after integration
        after = [self.system.storage.getParticle(i).radius for i in range(1,4)]

        # run checks (first particle excluded, hence it's radius should not change. The other should have changed, however, as they feel the thermostat on radius)
        self.assertEqual(before[0], after[0])
        self.assertNotEqual(before[1], after[1])
        self.assertNotEqual(before[2], after[2])

if __name__ == '__main__':
    unittest.main()
