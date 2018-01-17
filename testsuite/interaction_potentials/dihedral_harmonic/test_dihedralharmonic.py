#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import time
import math
import espressopp
import mpi4py.MPI as MPI

import unittest

def vec_pbc(u,v,cell):
   #vector u-v
   dx = u[0] - v[0]
   dy = u[1] - v[1]
   dz = u[2] - v[2]
   dx = dx - round(dx/cell[0])*cell[0]
   dy = dy - round(dy/cell[1])*cell[1] 
   dz = dz - round(dz/cell[2])*cell[2] 
   return dx,dy,dz

def abslen(u):
    return math.sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])

def cross(u,v):
    c=[0.0,0.0,0.0]
    c[0]=u[1]*v[2]-u[2]*v[1]
    c[1]=u[2]*v[0]-u[0]*v[2]
    c[2]=u[0]*v[1]-u[1]*v[0]
    return c

def dot(u,v):
    c=u[0]*v[0]+u[1]*v[1]+u[2]*v[2]
    return c

def calc_dihedral(self,quadrupleslist):
    #quadrupleslist contains [i,j,k,n]
    #returns torsional angle according to IUPAC convention (same convention as used in espressopp, dlpoly, gromacs, vmd, etc.)

    positions = []
    for pid in quadrupleslist:
      positions.append(self.system.storage.getParticle(pid).pos)
    
    rij = vec_pbc(positions[1],positions[0],self.box) 
    rjk = vec_pbc(positions[2],positions[1],self.box) 
    rkn = vec_pbc(positions[3],positions[2],self.box) 
    rijjk = cross(rij,rjk)
    rjkkn = cross(rjk,rkn)

    cos_phi = dot(rijjk,rjkkn)/(abslen(rijjk)*abslen(rjkkn))
    
    phi = math.acos(cos_phi)

    rcross = cross(rijjk,rjkkn) 
    signcheck = dot(rcross,rjk)
    if signcheck < 0.0: phi *= -1.0

    return phi

class TestDihedralHarmonic(unittest.TestCase):
    def setUp(self):

        system = espressopp.System()
        box = (10, 10, 10)
        self.box = box
        cutoff = 2.0
        skin = 1.0
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = skin
        system.comm = MPI.COMM_WORLD
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,cutoff,skin)
        cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, cutoff, skin)
        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
        self.system = system
 
    def test_phi0(self):
        #create system with three torsions, and for each one set a potential with a different phi0 values, to check that the interaction potential works for all combinations of positive and negative phi and phi0
       
        phi0 = [10.0,-170.0,170.0]

        particle_list = [
            #add particles with initial torsion angle +90
            (1, 0, espressopp.Real3D(2.0, 3.0, 3.0), 1.0),
            (2, 0, espressopp.Real3D(2.0, 2.0, 3.0), 1.0),
            (3, 0, espressopp.Real3D(2.0, 2.0, 2.0), 1.0),
            (4, 0, espressopp.Real3D(3.0, 2.0, 2.0), 1.0),
            #add particles with initial torsion angle 160
            (5, 0, espressopp.Real3D(2.0, 3.0, 3.0), 1.0),
            (6, 0, espressopp.Real3D(2.0, 2.0, 3.0), 1.0),
            (7, 0, espressopp.Real3D(2.0, 2.0, 2.0), 1.0),
            (8, 0, espressopp.Real3D(2.4, 0.8, 2.0), 1.0),
            #add particles with initial torsion angle -161
            ( 9, 0, espressopp.Real3D(2.0, 3.0, 3.0), 1.0),
            (10, 0, espressopp.Real3D(2.0, 2.0, 3.0), 1.0),
            (11, 0, espressopp.Real3D(2.0, 2.0, 2.0), 1.0),
            (12, 0, espressopp.Real3D(1.6, 0.8, 2.0), 1.0),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'pos', 'mass')
        self.system.storage.decompose()

        quadrupleslist = [[1,2,3,4],[5,6,7,8],[9,10,11,12]]
        torsiontuples = [(1,2,3,4),(5,6,7,8),(9,10,11,12)]
        bondtuples = [(1,2),(2,3),(3,4),(5,6),(6,7),(7,8),(9,10),(10,11),(11,12)]

        #add torsions
        interactions = []
        for i in xrange(3):
          fql = espressopp.FixedQuadrupleList(self.system.storage)
          fql.addQuadruples([torsiontuples[i]])
          interaction = espressopp.interaction.FixedQuadrupleListDihedralHarmonic(self.system,fql,potential=espressopp.interaction.DihedralHarmonic(K=1.0,phi0=phi0[i]*math.pi/180.0))
          self.system.addInteraction(interaction)
          interactions.append(interaction)

        #add bonds so that atoms in the torsions don't drift too far apart
        fpl = espressopp.FixedPairList(self.system.storage)
        fpl.addBonds(bondtuples)
        interaction = espressopp.interaction.FixedPairListHarmonic(self.system,fpl,potential=espressopp.interaction.Harmonic(K=1.0,r0=1.0))
        self.system.addInteraction(interaction)

        integrator = espressopp.integrator.VelocityVerlet(self.system)

        integrator.run(50)

        self.assertAlmostEqual(interactions[0].computeEnergy(),0.747885,places=5)
        self.assertAlmostEqual(interactions[1].computeEnergy(),0.099570,places=5)
        self.assertAlmostEqual(interactions[2].computeEnergy(),0.099570,places=5)

        self.assertAlmostEqual(calc_dihedral(self,quadrupleslist[0]),1.397549,places=5)
        self.assertAlmostEqual(calc_dihedral(self,quadrupleslist[1]),2.869874,places=5)
        self.assertAlmostEqual(calc_dihedral(self,quadrupleslist[2]),-2.869874,places=5)


if __name__ == '__main__':
    unittest.main()
