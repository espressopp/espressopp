#!/usr/bin/env python2
#  Copyright (C) 2016-2017(H)
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

###########################################################################
#  ESPResSo++                                                             #
#  Test script for Interpolation of tabulated LJ simulation               #
#                                                                         #
###########################################################################

import sys
import time
import os
import unittest
import mpi4py.MPI as MPI
import espressopp
import math
import numpy as np
import logging
from espressopp import Real3D, Int3D

# Input values for system
N = 10                                    # box size
size  = (float(N), float(N), float(N))
numParticles = 2                          # number of particles
nsnapshots = 20
nsteps  = 100                             # number of steps
skin    = 0.03                             # skin for Verlet lists
cutoff  = 0.2
splinetypes = ['Linear','Akima','Cubic']  # implemented spline types
tgt_probs = [0.55, 0.45]

######################################################################
## IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE   ##
######################################################################

# print '\n-- Tabulated Interpolation Test -- \n'
# print 'Steps: %3s' % nsteps
# print 'Particles: %3s' % numParticles
# print 'Cutoff: %3s' % cutoff

dim = 1
spline = 2
tabBonds = ["table_bi0.txt", "table_bi1.txt"]

# compute the number of cells on each node
def calcNumberCells(size, nodes, cutoff):
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

        system.storage.addParticle(0, Real3D(0.15, 0.1, 0))
        system.storage.addParticle(1, Real3D(0.4, 0.1, 0))

        system.storage.decompose()

        # integrator
        integrator = espressopp.integrator.VelocityVerlet(system)
        integrator.dt = 0.001

        # now build Verlet List
        # ATTENTION: you must not add the skin explicitly here
        logging.getLogger("Interpolation").setLevel(logging.INFO)
        vl = espressopp.VerletList(system, cutoff = cutoff)

        scaling_cv = 1.5
        alpha      = 0.1

        # ATTENTION: auto shift was enabled
        # bonds
        writeTabFile(tabBonds[0], k=100000, x0=0.26, N=500, low=0.01, high=0.6)
        writeTabFile(tabBonds[1], k= 50000, x0=0.24, N=500, low=0.01, high=0.6)

        fpl = espressopp.FixedPairList(system.storage)
        fpl.addBonds([(0,1)])
        potBond = espressopp.interaction.TabulatedSubEns()
        potBond.addInteraction(1, tabBonds[0],
            espressopp.RealND([0.262, 0., 0., scaling_cv*0.424, 0., 0.]))
        potBond.addInteraction(1, tabBonds[1],
            espressopp.RealND([0., 0., 0., 0., 0., 0.]))
        potBond.alpha_set(alpha)
        cv_bl = espressopp.FixedPairList(system.storage)
        cv_bl.addBonds([])
        potBond.colVarBondList = cv_bl
        cv_al = espressopp.FixedTripleList(system.storage)
        cv_al.addTriples([])
        potBond.colVarAngleList = cv_al

        # Renormalize the CVs
        potBond.colVarSd_set(0,   0.0119)

        # Target probabilities
        potBond.targetProb_set(0, tgt_probs[0])
        potBond.targetProb_set(1, tgt_probs[1])

        interBond = espressopp.interaction.FixedPairListTabulatedSubEns(system, fpl, potBond)
        system.addInteraction(interBond)

        temp = espressopp.analysis.Temperature(system)

        temperature = temp.compute()
        Ek = 0.5 * temperature * (2 * numParticles)
        Ep = interBond.computeEnergy()

        # langevin thermostat
        langevin = espressopp.integrator.LangevinThermostat(system)
        integrator.addExtension(langevin)
        langevin.gamma = 1.0
        langevin.temperature = 2.479 # in kJ/mol

        # sock = espressopp.tools.vmd.connect(system)

        configurations = espressopp.analysis.Configurations(system)
        self.configuration = configurations
        self.system = system
        self.temp = temp
        self.integrator = integrator
        self.interBond = interBond
        self.potBond = potBond


class TestSubEnsBond(makeConf):
    def test_sub_ens_bond(self):
        weights = []
        for k in range(nsnapshots):
            self.integrator.run(nsteps)
            # espressopp.tools.vmd.imd_positions(system, sock)
            temperature = self.temp.compute()
            Ek = 0.5 * temperature * (2 * numParticles)
            Eb = self.interBond.computeEnergy()
            wts = self.potBond.weight_get()
            cv = self.potBond.colVar
            weights.append([wts[0],wts[1]])
            if k % 10 == 0:
                print 'Step %6d:' % ((k+10)*nsteps),
                print "Weights: (%3f, %3f)" % (wts[0], wts[1]),
                print "avg: %3f" % np.mean([w[0] for w in weights]),
                print "ColVar: (%3f)" % (cv[0])
        self.assertNotEqual(np.mean([w[0] for w in weights]), 0)
        self.assertNotEqual(np.mean([w[1] for w in weights]), 0)

if __name__ == '__main__':
    unittest.main()
