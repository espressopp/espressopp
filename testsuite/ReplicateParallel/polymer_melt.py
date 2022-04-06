#!/usr/bin/env python3
#  Copyright (C) 2012-2017(H)
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
#                                                                         #
#  ESPResSo++ Python script for a Polymer Melt System including           #
#  runtime details
#                                                                         #
###########################################################################

import time
import espressopp
import unittest
import numpy as np

nsteps      = 10
isteps      = 100
rc          = pow(2.0, 1.0/6.0)
skin        = 0.4
timestep    = 0.005

# set temperature to None for NVE-simulations
temperature = 1.0

repl = (2,2,2)

######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################
print(espressopp.Version().info())
print('Setting up simulation ...')

def init_polymer_melt_simulation(use_replicate_parallel):

    start_time = time.process_time()

    bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.lammps.read('polymer_melt.lammps')

    if use_replicate_parallel:
        # rp stores bonds, angles, and positions
        rp = espressopp.tools.ReplicateParallel()
        num_particles, Lx, Ly, Lz = rp.replicate(bonds, angles, x, y, z, Lx, Ly, Lz, *repl)
    else:
        bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.replicate(bonds, angles, x, y, z, Lx, Ly, Lz, *repl)
        num_particles = len(x)

    density = num_particles / (Lx * Ly * Lz)
    box = (Lx, Ly, Lz)
    system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)

    if use_replicate_parallel:
        props = ['type', 'mass']
        num_particles_seed = len(x)
        seed_particles = []
        for i in range(num_particles_seed):
            part = [0, 1.0]
            seed_particles.append(part)
        # rp workers add replicated particles in parallel
        rp.addParticles(system.storage, 1, seed_particles, *props)
    else:
        # add particles to the system and then decompose
        # do this in chunks of 1000 particles to speed it up
        props = ['id', 'type', 'mass', 'pos']
        new_particles = []
        for i in range(num_particles):
            part = [i + 1, 0, 1.0, espressopp.Real3D(x[i], y[i], z[i])]
            new_particles.append(part)
            if i % 1000 == 0:
                system.storage.addParticles(new_particles, *props)
                system.storage.decompose()
                new_particles = []
        system.storage.addParticles(new_particles, *props)

    system.storage.decompose()

    # Lennard-Jones with Verlet list
    vl      = espressopp.VerletList(system, cutoff = rc)
    potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
    interLJ = espressopp.interaction.VerletListLennardJones(vl)
    interLJ.setPotential(type1=0, type2=0, potential=potLJ)
    system.addInteraction(interLJ)

    # FENE bonds
    fpl = espressopp.FixedPairList(system.storage)
    if use_replicate_parallel:
        # rp workers add replicated bonds in parallel
        rp.addBonds(fpl)
    else:
        fpl.addBonds(bonds)
    potFENE = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
    interFENE = espressopp.interaction.FixedPairListFENE(system, fpl, potFENE)
    system.addInteraction(interFENE)

    # Cosine with FixedTriple list
    ftl = espressopp.FixedTripleList(system.storage)
    if use_replicate_parallel:
        # rp workers add replicated triples in parallel
        rp.addTriples(ftl)
    else:
        ftl.addTriples(angles)
    potCosine = espressopp.interaction.Cosine(K=1.5, theta0=3.1415926)
    interCosine = espressopp.interaction.FixedTripleListCosine(system, ftl, potCosine)
    system.addInteraction(interCosine)

    end_time = time.process_time()

    # get positions, bonds and angles
    configurations = espressopp.analysis.Configurations(system, pos=True, vel=False, force=False)
    configurations.gather()
    pos = [configurations[0][i] for i in range(num_particles)]
    bonds = np.array(fpl.getBonds())
    angles = np.array(ftl.getTriples())

    return end_time - start_time, pos, bonds, angles

class TestReplicateParallel(unittest.TestCase):

    def test1(self):
        ''' Ensure that particles, bonds and triples are the same '''

        t0, pos0, bonds0, angles0 = init_polymer_melt_simulation(False)
        t1, pos1, bonds1, angles1 = init_polymer_melt_simulation(True)

        print(f"Initialize polymer melt simulation using replicate:         {t0:.2f} seconds")
        print(f"Initialize polymer melt simulation using ReplicateParallel: {t1:.2f} seconds")

        self.assertEqual(len(pos0), len(pos1))
        diff = [(pos0[i]-pos1[i]).sqr() for i in range(len(pos1))]
        for d in diff:
            self.assertAlmostEqual(d,0.0,8)

        np.testing.assert_equal(bonds0, bonds1)

        np.testing.assert_equal(angles0, angles1)

if __name__ == "__main__":
    unittest.main()