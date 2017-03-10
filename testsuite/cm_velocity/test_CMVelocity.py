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
