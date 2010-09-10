from espresso import *
from espresso import unittest
import MPI

class TestParticleLocal(unittest.TestCase) :
    def test0get(self):
        system = System()
        system.rng = esutil.RNG()
        system.bc = bc.OrthorhombicBC(system.rng, (10.0, 10.0, 10.0))
        system.storage = espresso.storage.DomainDecomposition(
            system=system, comm=MPI.COMM_WORLD, 
            nodeGrid=(1,1,1), cellGrid=(1,1,1))
        p = system.storage.addParticle(0, (1.0, 1.0, 1.0))
        p.v = Real3D(1.0, 1.0, 1.0)
        self.assertAlmostEqualReal3D(p.v, Real3D(1.0, 1.0, 1.0))

if __name__ == "__main__":
    unittest.main()
