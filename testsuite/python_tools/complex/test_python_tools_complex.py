#!/usr/bin/env python
# -*- coding: utf-8 -*-

#unittests for python tools in $ESPRESSOPPHOME/tools
#this file uses a system containing several particles, for testing more complex tools 
#input/output tools should be tested using test_python_tools_io.py


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
            (1, 1,  0.1, espressopp.Real3D(3.0, 1.0, 4.0), 2.0),
            (2, 1, -0.5, espressopp.Real3D(2.0, 2.0, 4.0), 2.0),
            (3, 1, -0.5, espressopp.Real3D(1.0, 1.0, 4.5), 2.0),
            (4, 1,  0.5, espressopp.Real3D(4.0, 1.0, 4.0), 2.0),
            (5, 1, -0.5, espressopp.Real3D(5.0, 2.0, 4.0), 2.0),
            (6, 1,  0.2, espressopp.Real3D(6.0, 1.0, 4.5), 2.0),
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'q', 'pos', 'mass')
        self.system.storage.decompose()

    def test_boresch_restraints(self):

        restraintAtoms = {'A':1,'B':2,'C':3,'a':4,'b':5,'c':6} #A,B,C in ligand, a,b,c in protein, original atomistic indices
        restraintK = {'aA':4184,'baA':41.84,'aAB':41.84,'aABC':41.84,'cbaA':41.84,'baAB':41.84} #kJ mol-1 nm-2,kJ mol-1 rad-2, all 10 kcal as in JPCB 2003
        restraintR0 = {'aA':1.01,'baA':120.0,'aAB':120.0,'aABC':100.0,'cbaA':-170.0,'baAB':-105.0} #nm, degrees

        restraint_interactions = espressopp.tools.applyBoreschRestraints(self.system,restraintAtoms,restraintK,restraintR0)
        dhdlRstr = 0.0
        for rt in restraint_interactions.values(): dhdlRstr+=rt.computeEnergy() #energyDeriv = energy for restraints
        self.assertAlmostEqual(dhdlRstr,141.92348724,places=5)

    def test_self_excl_energy(self):

        exclusions = []
        energy = espressopp.tools.energy.getSelfExclEnergyReactionField(self.system,exclusions,prefactor=1.0,epsilon1=1,epsilon2=80,rc=3.0,pidlist=[3,4,5,6])
        self.assertAlmostEqual(energy,-0.1962732919254,places=5)
        energy = espressopp.tools.energy.getSelfExclEnergyReactionField(self.system,exclusions,prefactor=1.0,epsilon1=1,epsilon2=80,rc=3.0,nParticles=6)
        self.assertAlmostEqual(energy,-0.2608695652173,places=5)


if __name__ == '__main__':
    unittest.main()
