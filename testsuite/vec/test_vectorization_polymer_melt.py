#!/usr/bin/env python3
#  Copyright (C) 2021
#      Max Planck Institute for Polymer Research & JGU Mainz
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

import unittest
import espressopp
from espressopp.tools import readxyz
import time

def generate_md(use_vec=True):
    print('{}USING VECTORIZATION'.format('NOT ' if not use_vec else ''))
    nsteps      = 1
    isteps      = 10
    repl        = (1,1,1)
    #
    # NOTE: For performance comparison increase isteps to 1000 and repl to (2,2,2)
    #
    rc          = pow(2.0, 1.0/6.0)
    skin        = 0.4
    timestep    = 0.005

    # ensure deterministic trajectories
    temperature = None

    bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.lammps.read('polymer_melt.lammps')
    bonds, angles, x, y, z, Lx, Ly, Lz = espressopp.tools.replicate(bonds, angles, x, y, z, Lx, Ly, Lz, xdim=repl[0], ydim=repl[1], zdim=repl[2])
    box = (Lx, Ly, Lz)
    num_particles = len(x)
    if use_vec:
        system, integrator = espressopp.vec.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)
    else:
        system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)

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

    if use_vec:
        # Lennard-Jones with Verlet list
        vl      = espressopp.vec.VerletList(system, cutoff = rc)
        potLJ   = espressopp.vec.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
        interLJ = espressopp.vec.interaction.VerletListLennardJones(vl)
        interLJ.setPotential(type1=0, type2=0, potential=potLJ)
        system.addInteraction(interLJ)

        # FENE bonds
        fpl       = espressopp.vec.FixedPairList(system.storage)
        fpl.addBonds(bonds)
        potFENE   = espressopp.vec.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
        interFENE = espressopp.vec.interaction.FixedPairListFENE(system, fpl, potFENE)
        system.addInteraction(interFENE)

        # Cosine with FixedTriple list
        ftl         = espressopp.vec.FixedTripleList(system.storage)
        ftl.addTriples(angles)
        potCosine   = espressopp.vec.interaction.Cosine(K=1.5, theta0=3.1415926)
        interCosine = espressopp.vec.interaction.FixedTripleListCosine(system, ftl, potCosine)
        system.addInteraction(interCosine)
    else:
        # Lennard-Jones with Verlet list
        vl      = espressopp.VerletList(system, cutoff = rc, useBuffers=False)
        potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
        interLJ = espressopp.interaction.VerletListLennardJones(vl)
        interLJ.setPotential(type1=0, type2=0, potential=potLJ)
        system.addInteraction(interLJ)

        # FENE bonds
        fpl       = espressopp.FixedPairList(system.storage)
        fpl.addBonds(bonds)
        potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
        interFENE = espressopp.interaction.FixedPairListFENE(system, fpl, potFENE)
        system.addInteraction(interFENE)

        # Cosine with FixedTriple list
        ftl         = espressopp.FixedTripleList(system.storage)
        ftl.addTriples(angles)
        potCosine   = espressopp.interaction.Cosine(K=1.5, theta0=3.1415926)
        interCosine = espressopp.interaction.FixedTripleListCosine(system, ftl, potCosine)
        system.addInteraction(interCosine)

    print('')
    print('number of particles = ', num_particles)
    print("storage             = ", system.storage.__class__.__name__)
    print("integrator          = ", integrator.__class__.__name__)
    print("verletlist          = ", ".".join([vl.__class__.__module__,vl.__class__.__name__]))
    print("interaction         = ", ".".join([interLJ.__class__.__module__,interLJ.__class__.__name__]))
    print('')

    espressopp.tools.analyse.info(system, integrator)

    start_time = time.process_time()
    for k in range(nsteps):
        integrator.run(isteps)
        espressopp.tools.analyse.info(system, integrator)
    end_time = time.process_time()

    espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)

    # retrieve particle positions after run
    configurations = espressopp.analysis.Configurations(system, pos=True, vel=True, force=True)
    configurations.gather()

    return [configurations[0][i] for i in range(num_particles)]

class TestVectorization(unittest.TestCase):

    def test1(self):
        ''' Ensure that positions after integration are the same for both vec and non-vec versions '''
        print('-'*70)
        pos0 = generate_md(True)
        print('-'*70)
        pos1 = generate_md(False)
        print('-'*70)

        self.assertEqual(len(pos0), len(pos1))
        diff = [(pos0[i]-pos1[i]).sqr() for i in range(len(pos1))]
        for d in diff:
            self.assertAlmostEqual(d,0.0,8)

if __name__ == "__main__":
    unittest.main()
