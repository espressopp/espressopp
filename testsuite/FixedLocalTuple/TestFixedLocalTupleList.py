#!/usr/bin/env python2
#
#  Copyright (C) 2013-2017(H)
#      Max Planck Institute for Polymer Research
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
# 
# -*- coding: utf-8 -*-
#

import espressopp
import unittest
import mpi4py.MPI as MPI
import math

from espressopp import Real3D

class TestFixedLocalTupleList(unittest.TestCase) :

    def setUp(self) :
        system = espressopp.System()
        
        rng  = espressopp.esutil.RNG()
        
        N    = 4
        SIZE = float(N)
        box  = Real3D(SIZE)
        bc   = espressopp.bc.OrthorhombicBC(None, box)
        
        system.bc = bc
        
        # a small skin avoids rounding problems
        
        system.skin = 0.001
    
        cutoff = SIZE/2. - system.skin 
        
        comm = espressopp.MPI.COMM_WORLD
        
        nodeGrid       = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,cutoff,system.skin)
        cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, cutoff, system.skin)
        
        print 'NodeGrid = %s'%(nodeGrid,)
        print 'CellGrid = %s'%cellGrid

        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
        pid = 0

        for i in xrange(N):
            for j in xrange(N):
                for k in xrange(N):
                    
                    r = 0.5
                    x = (i + r) / N * SIZE
                    y = (j + r) / N * SIZE
                    z = (k + r) / N * SIZE
   
                    system.storage.addParticle(pid, Real3D(x, y, z))

                    pid = pid + 1

        for i in xrange(N):
            for j in xrange(N):
                for k in xrange(N):
                    
                    r = 0.25
                    x = (i + r) / N * SIZE
                    y = (j + r) / N * SIZE
                    z = (k + r) / N * SIZE
   
                    system.storage.addParticle(pid, Real3D(x, y, z))

                    pid = pid + 1
                    
        system.storage.decompose()
        
	# now build Fixed Local Tuple List
        tuplelist = espressopp.FixedLocalTupleList(system.storage)

	self.system = system
	self.N = N
	self.tuplelist = tuplelist

    # This function checks the size of empty FixedLocalTupleList.
    def test_create_fixedtuplelist(self) :
        self.assertEqual(sum(self.tuplelist.size(), 0), 0)

    # This function checks the python interface of FixedLocalTupleList.
    # For addTuple() and size(), this function test
    # whether the added tuple number equals the return value of size().
    # For getTuples(), this function test
    # whether a tuplelist obtained by getTuples equals added tuplelist
    def test_add_get_fixedtuplelist(self) :
	system = self.system
	N = self.N
	tuplelist = self.tuplelist

        # FixedLocalTupleList contain particles
        num_constrain = N*N
        stored = []
        for i in range(N*N*N/num_constrain):
            tuple = []
            for j in range(num_constrain):
                tuple.append(num_constrain*i + j)
            tuplelist.addTuple(tuple)
            stored.append(tuple)

        num_constrain = N*N/2
        for i in range(N*N*N/num_constrain, 2*N*N*N/num_constrain):
            tuple = []
            for j in range(num_constrain):
                tuple.append(num_constrain*i + j)
            tuplelist.addTuple(tuple)
            stored.append(tuple)

	# check the size of FixedLocalTupleList
        self.assertEqual(sum(tuplelist.size(), 0), 1.5*N*N*N/num_constrain)

        # check the contained particles id
        g_tuplelist = tuplelist.getTuples()
        s_id = 0
        for i in range(3*N*N*N/num_constrain/2):
            for j in range(espressopp.MPI.COMM_WORLD.size):
                if stored[s_id] in g_tuplelist[j]:
                    break
            self.assertEqual(stored[s_id] in g_tuplelist[j], True)
            s_id += 1

if __name__ == "__main__":
    unittest.main()

