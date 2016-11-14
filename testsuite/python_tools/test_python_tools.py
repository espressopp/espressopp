#!/usr/bin/env python
# -*- coding: utf-8 -*-

#unittests for python tools in $ESPRESSOPPHOME/tools

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
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, 1.5, 0.3)
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
        nparticles = 2
        mass,charge,index,name,ptype = espressopp.tools.readSimpleSystem(filename,nparticles,header=1)

        self.assertEqual(mass[1],5.0)
        self.assertEqual(charge[1],-0.7)
        self.assertEqual(index[1],8)
        self.assertEqual(name[1],'name2')
        self.assertEqual(ptype[1],'type2')


if __name__ == '__main__':
    unittest.main()
