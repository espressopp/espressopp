#!/usr/bin/env python2

import sys
import unittest
import time
import os
import mpi4py.MPI as MPI
import espressopp
import math
import numpy as np
import logging
from espressopp import Real3D, Int3D

print "Comparison Tabulated dihedral vs Harmonic dihedral"

# Input values for system
N = 10                                    # box size
size  = (float(N), float(N), float(N))
numParticles = 4                          # number of particles
nsnapshots = 100
nsteps  = 10                                # number of steps
skin    = 0.03                             # skin for Verlet lists
cutoff  = 0.2
splinetypes = ['Linear','Akima','Cubic']  # implemented spline types

dim = 1
spline = 2

tabBond = "table_b0.txt"
tabDih  = "table_d_h10.txt"

print "Max simulation time: {:d}".format(nsnapshots * nsteps)

def calcNumberCells(size, nodes, cutoff):
    # compute the number of cells on each node
    ncells = 1
    while size / (ncells * nodes) >= cutoff:
       ncells = ncells + 1
    return ncells - 1

# writes the tabulated file
def writeTabFile(name, k, x0, N, low=0.0, high=2.5):
    outfile = open(name, "w")
    delta = (high - low) / (N - 1)

    for i in range(N):
        r = low + i * delta
        energy = 0.5*k*(r-x0)**2
        force  = -k*(r-x0)
        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))

    outfile.close()

class makeConf(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        system.rng  = espressopp.esutil.RNG()
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
        system.skin = skin

        comm = MPI.COMM_WORLD

        nodeGrid = Int3D(1, 1, comm.size)
        cellGrid = Int3D(
            calcNumberCells(size[0], nodeGrid[0], cutoff),
            calcNumberCells(size[1], nodeGrid[1], cutoff),
            calcNumberCells(size[2], nodeGrid[2], cutoff)
            )
        system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

        pid = 0

        system.storage.addParticle(0, Real3D(6.1162456968, 3.9088541374, 9.4409851324))
        system.storage.addParticle(1, Real3D(5.9516473019, 4.2535440563, 9.4778266528))
        system.storage.addParticle(2, Real3D(5.8160837525, 4.1043280354, 9.8101428319))
        system.storage.addParticle(3, Real3D(6.0145256344, 4.0133146160, 9.4982604364))

        system.storage.decompose()
        self.system = system

class TestTabDih(makeConf):
    def test_tab_dih(self):
        # integrator
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.001

        # now build Verlet List
        # ATTENTION: you must not add the skin explicitly here
        logging.getLogger("Interpolation").setLevel(logging.INFO)
        vl = espressopp.VerletList(self.system, cutoff = cutoff)


        # bonds
        writeTabFile(tabBond, k=400000, x0=0.4, N=500, low=0.01, high=0.80)

        fpl = espressopp.FixedPairList(self.system.storage)
        fpl.addBonds([(0,1),(1,2),(2,3)])
        potBond = espressopp.interaction.Tabulated(itype=spline, filename=tabBond)
        interBond = espressopp.interaction.FixedPairListTabulated(self.system, fpl, potBond)
        self.system.addInteraction(interBond)

        # dihedral
        spring = 10
        x_rest = 0.0
        writeTabFile(tabDih, k=spring, x0=x_rest, N=500, low=-np.pi, high=np.pi)

        fql = espressopp.FixedQuadrupleList(self.system.storage)
        fql.addQuadruples([(0,1,2,3)])
        potDihed1 = espressopp.interaction.TabulatedDihedral(itype=spline,
                                            filename=tabDih)
        interDihed1 = espressopp.interaction.FixedQuadrupleListTabulatedDihedral(self.system, fql, potDihed1)
        potDihed2 = espressopp.interaction.DihedralHarmonic(spring, x_rest)
        interDihed2 = espressopp.interaction.FixedQuadrupleListDihedralHarmonic(self.system, fql, potDihed2)
        self.system.addInteraction(interDihed1)


        temp = espressopp.analysis.Temperature(self.system)

        temperature = temp.compute()
        Ek = 0.5 * temperature * (3 * numParticles)
        Ep = interDihed1.computeEnergy()

        # langevin thermostat
        langevin = espressopp.integrator.LangevinThermostat(self.system)
        integrator.addExtension(langevin)
        langevin.gamma = 1.0
        langevin.temperature = 2.479 # in kJ/mol
        print "Running at temperature T = {:.3f} kJ/mol/k_B".format(langevin.temperature)


        start_time = time.clock()
        print " ***"
        print "{:8s} {:8s} {:8s}".format("Step","E_tab","E_harm")
        for k in range(nsnapshots):
            Ed1, Ed2 = interDihed1.computeEnergy(), interDihed2.computeEnergy()
            if k % 10 == 0:
                self.assertAlmostEqual(Ed1, Ed2, places=2)
                print '{:8d} {:8f} {:8f}'.format(((k+10)*nsteps),Ed1, Ed2)
            integrator.run(nsteps)


if __name__ == '__main__':
    unittest.main()
