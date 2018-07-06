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
#unittests for python tools in $ESPRESSOPPHOME/tools
#this file uses a system containing only two particles, for testing input/output tools with minimum computational time

import sys
import time
import espressopp
import mpi4py.MPI as MPI
import unittest


class TestPythonTools(unittest.TestCase):

    def setUp(self):
        system = espressopp.System()
        box = (10, 10, 10)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc=1.5,skin=system.skin)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid,rc= 1.5, skin=system.skin)
        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
        self.system = system

        #add particles
        particle_list = [
            (1, 1,  0.178, espressopp.Real3D(2.0, 3.0, 4.0), 1.008),
            (2, 2, -0.513, espressopp.Real3D(2.0, 3.0, 4.0), 41.54),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass')
        self.system.storage.decompose()

    def test_psf_read_write(self):

        filename = 'temp.psf'
        typenames = {1: 'AA', 2: 'BB'}

        espressopp.tools.psfwrite(filename, self.system, maxdist=None, molsize=0, typenames=typenames)
        pid,segname,resindex,resname,atomname,atomtype,mass,charge = espressopp.tools.psfread(filename)

        self.assertEqual(pid[0],1)
        self.assertEqual(pid[1],2)
        self.assertEqual(atomtype[0].strip(),'AA')
        self.assertEqual(atomtype[1].strip(),'BB')
        self.assertEqual(mass[0],1.008)
        self.assertEqual(mass[1],41.54)
        self.assertEqual(charge[0],0.178)
        self.assertEqual(charge[1],-0.513)

    def test_pdb_read_write(self):

        filename = 'temp.pdb'
        typenames = {1: 'AA', 2: 'BB'}

        espressopp.tools.pdbwrite(filename, self.system, molsize=0, append=False, typenames=typenames)
        index,atomname,resname,resid,x,y,z,alpha,beta,segid,element = espressopp.tools.pdbread(filename,natoms=2,header=2)

        self.assertEqual(index[0],1)
        self.assertEqual(index[1],2)
        self.assertEqual(atomname[0].strip(),'AA')
        self.assertEqual(atomname[1].strip(),'BB')
        self.assertEqual(resid[0],0)
        self.assertEqual(resid[1],0)
        self.assertEqual(x[0],2.0)
        self.assertEqual(y[0],3.0)
        self.assertEqual(z[0],4.0)
        self.assertEqual(x[1],2.0)
        self.assertEqual(y[1],3.0)
        self.assertEqual(z[1],4.0)

    def test_read_simple_system(self):

        filename = 'test.top'
        ff = open(filename,'w')
        ff.write('#mass,charge,index,name,ptype\n')
        ff.write('1.0 -1.0 4 name1 type1\n')
        ff.write('5.0 -0.7 8 name2 type2\n')
        ff.close()

        nparticles = 2
        mass,charge,index,name,ptype = espressopp.tools.readSimpleSystem(filename,nparticles,header=1)

        self.assertEqual(mass[1],5.0)
        self.assertEqual(charge[1],-0.7)
        self.assertEqual(index[1],8)
        self.assertEqual(name[1],'name2')
        self.assertEqual(ptype[1],'type2')


if __name__ == '__main__':
    unittest.main()
