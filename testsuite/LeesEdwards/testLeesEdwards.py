#!/usr/bin/env python3
#  Copyright (C) 2022
#      Data Center, Johannes Gutenberg University Mainz
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

import sys
import mpi4py.MPI as MPI
import espressopp
from espressopp import Real3D, Int3D
import time
import math
import numpy as np

import unittest

class TestLeesEdwards(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        rng = espressopp.esutil.RNG()
        rng.seed(1)
        system.rng = rng
        box = (10, 10, 10)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.5
        system.comm = MPI.COMM_WORLD
        self.system = system
    
    def test_normal(self):
        # set up normal domain decomposition
        box=(10,10,10)
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box, rc=1.5, skin=0.3)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc=1.5, skin=0.3)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)
        
        # add some particles (normal, coarse-grained particles only)
        x = [ espressopp.Real3D(5.0, 5.0, 5.0),
              espressopp.Real3D(5.0, 5.0, 5.0) ]
        
        v = [ espressopp.Real3D(0.0,0.0,1.0),
              espressopp.Real3D(0.0,0.0,-1.0) ]
        
        particle_list = [
            (1, 0, x[0], v[0], 1.0),
            (2, 0, x[1], v[1], 1.0)
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'pos', 'v', 'mass')
        self.system.storage.decompose()
        
        # generate a verlet list
        vl = espressopp.VerletList(self.system, cutoff=1.5)
        
        # integrator
        integrator = espressopp.integrator.VelocityVerletLE(self.system,shear=0.75)
        integrator.dt = 0.001
        
        for i in range(100):
            integrator.run(50)
            #espressopp.tools.xyzfilewrite("out.xyz", self.system, velocities = False, charge = False, append=True, atomtypes={0:'X'})
            #espressopp.tools.pdbwrite("out.pdb", self.system, molsize=2,append=True)
        
        #conf  = espressopp.analysis.Configurations(self.system)
        #conf.capacity=2
        #conf.gather()
        p1=self.system.storage.getParticle(1).pos
        p2=self.system.storage.getParticle(2).pos
        
        #print("RES: ",p1,"|",p2)
        #print("VEL: ",self.system.storage.getParticle(1).v,"|",self.system.storage.getParticle(2).v)
        
        
        #the checks follow
        
        #Check if the forces are (almost) equal
        #self.assertEqual(p1[0],p2[0])
        self.assertAlmostEqual(p1[0],p2[0],places=2)

if __name__ == '__main__':
    unittest.main()
