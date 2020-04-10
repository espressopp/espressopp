#!/usr/bin/env python2
#
#  Copyright (C) 2020(Dr. Horacio Andres Vargas Guzman)
#      Institute Jozef Stefan
#      Max Planck Institute for Polymer Research
#
#  This file is part of ESPResSo++ -> Powered by HeSpaDDA algorithm developed by horacio.v.g@gmail.com
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
#
# -*- coding: utf-8 -*-

import espressopp
import mpi4py.MPI as MPI
from espressopp import Int3D
from espressopp import Real3D

import unittest
import os
import shutil

box_size=[40.0, 20.0, 20.0]
rc=1.12
skin=0.3
pTotal=[16, 64, 128, 250, 256, 384, 512, 1000, 1024, 2000]
pX=[0,0,0,0,0,0,0,0,0,0]
pY=[0,0,0,0,0,0,0,0,0,0]
pZ=[0,0,0,0,0,0,0,0,0,0]

# expected outputs
hdd1=[4,2,2]
hdd2=[16,2,2]
hdd3=[8,4,4]
hdd4=[10,5,5]
hdd5=[16,4,4]
hdd6=[24,4,4]
hdd7=[32,4,4]
hdd8=[40,5,5]
hdd9=[16,8,8]
hdd10=[20,10,10]

# AdResS setup
eh_size=10.99  # ex_size+hy_size
ratioMS=3.0
idealGas=0   # No load at all in the low-res region
slabMS=[1,0,0]

class TestOneScale(unittest.TestCase):
	def test_OneScaleHDD(self):
		print "Checking the domain decomposition for a heterogeneous (oneScale) MD simulations"
		# set bassic parameters
		global box_size,rc,skin
		for i in range(1,11):
			pX[i-1], pY[i-1], pZ[i-1] = espressopp.tools.decomp.nodeGrid(pTotal[i-1],box_size,rc,skin)
		self.assertAlmostEqual(pX[0],hdd1[0],places=2)
		self.assertAlmostEqual(pY[9],hdd10[1],places=2)
		self.assertAlmostEqual(pZ[5], hdd6[2], places=2)

if __name__ == '__main__':
    unittest.main()		
