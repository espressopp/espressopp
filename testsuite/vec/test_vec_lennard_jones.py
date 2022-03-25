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
    #
    # NOTE: For performance comparison increase isteps to 1000
    #
    rc          = 2.5
    skin        = 0.3
    timestep    = 0.005
    epsilon     = 1.0
    sigma       = 1.0

    # ensure deterministic trajectories
    temperature = None

    xyz_file = "lennard_jones_fluid_10000_2048.xyz"
    pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = readxyz(xyz_file)
    box = (Lx, Ly, Lz)
    num_particles = len(pid)

    if use_vec:
        system, integrator = espressopp.vec.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)
    else:
        system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)

    props = ['id', 'type', 'mass', 'pos', 'v']
    new_particles = []
    for i in range(num_particles):
        part = [i + 1, 0, 1.0, espressopp.Real3D(x[i], y[i], z[i]), espressopp.Real3D(vx[i], vy[i], vz[i])]
        new_particles.append(part)
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()

    # Lennard-Jones with Verlet list
    if use_vec:
        vl      = espressopp.vec.VerletList(system, cutoff = rc)
        interLJ = espressopp.vec.interaction.VerletListLennardJones(vl)
        potLJ   = espressopp.vec.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)
    else:
        vl      = espressopp.VerletList(system, cutoff = rc)
        interLJ = espressopp.interaction.VerletListLennardJones(vl)
        potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)

    interLJ.setPotential(type1=0, type2=0, potential=potLJ)
    system.addInteraction(interLJ)

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
