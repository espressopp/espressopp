import espressopp
import math
import unittest

# initial parameters of the simulation
L              = 5
box            = (L, L, L)
rc             = 1.5

sigma          = 1.
phi            = 0.

class makeConf(unittest.TestCase):
    def setUp(self):

        system, integrator = espressopp.standard_system.Default(box, rc=rc, skin=0.3, dt=0.005, temperature=1.)
        system.storage.addParticles([[1,0,espressopp.Real3D(0,0,1)]],'id','type','pos')
        system.storage.addParticles([[2,0,espressopp.Real3D(0,0,1+math.pow(2,1./6.))]],'id','type','pos')
        system.storage.decompose()

        # non-bonded LJcos potential
        vl = espressopp.VerletList(system, cutoff=rc)
        LJcos = espressopp.interaction.LJcos(phi=phi)
        LJcosInter = espressopp.interaction.VerletListLJcos(vl)
        LJcosInter.setPotential(type1=0, type2=0, potential=LJcos)
        system.addInteraction(LJcosInter)

        # set self
        self.system = system
        self.LJcosInter = LJcosInter

class TestLJcos(makeConf):
    def test_ljcos(self):

        for phi in range (0, 11):
            phi = 0.1 * phi
            # update potential
            LJcos = espressopp.interaction.LJcos(phi=phi)
            self.LJcosInter.setPotential(type1=0, type2=0, potential=LJcos)

            Epot = espressopp.analysis.EnergyPot(self.system)
            print 'phi =', phi, 'Epot = ', Epot.compute()

            self.assertAlmostEqual(Epot.compute(), -phi, places=10)

if __name__ == '__main__':
    unittest.main()

