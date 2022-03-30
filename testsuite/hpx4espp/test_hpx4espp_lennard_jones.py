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
from espressopp import hpx4espp
from espressopp.tools import readxyz
import time
import os, sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--xyz_file", type=str, default="lennard_jones_fluid_10000_2048.xyz")
args, _ = parser.parse_known_args()
xyz_file = args.xyz_file
assert os.path.isfile(args.xyz_file), f"xyz_file={args.xyz_file} does not exist. " \
    "Provide the path to file using the --xyz_file [XYZ_FILE] argument. The file is " \
    "testsuite/vec/lennard_jones_fluid_10000_2048.xyz in the source directory."

def generate_md(use_hpx=True):
    print('{}USING HPX4ESPP'.format('NOT ' if not use_hpx else ''))

    if use_hpx:
        Default                      = hpx4espp.standard_system.Default
        VerletList                   = hpx4espp.VerletList
        VerletListLennardJones       = hpx4espp.interaction.VerletListLennardJones
        LennardJones                 = hpx4espp.interaction.LennardJones
    else:
        Default                      = espressopp.standard_system.Default
        VerletList                   = espressopp.VerletList
        VerletListLennardJones       = espressopp.interaction.VerletListLennardJones
        LennardJones                 = espressopp.interaction.LennardJones

    # nsteps      = 1
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

    pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = readxyz(args.xyz_file)
    box = (Lx, Ly, Lz)
    num_particles = len(pid)

    numSubs = {"numSubs":8} if use_hpx else {}

    system, integrator = Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature, **numSubs)

    props = ['id', 'type', 'mass', 'pos', 'v']
    new_particles = []
    for i in range(num_particles):
        part = [i + 1, 0, 1.0, espressopp.Real3D(x[i], y[i], z[i]), espressopp.Real3D(vx[i], vy[i], vz[i])]
        new_particles.append(part)
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()

    # Lennard-Jones with Verlet list
    vl      = VerletList(system, cutoff = rc)
    interLJ = VerletListLennardJones(vl)
    potLJ   = LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0)

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

    hpx = hpx4espp.HPXRuntime()

    if use_hpx:
        hpx.start(threads=4)
    start_time = time.process_time()

    integrator.run(isteps)

    end_time = time.process_time()
    if use_hpx:
        hpx.stop()

    espressopp.tools.analyse.info(system, integrator)
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
    # do not pass cmd-line arguments to unittest
    sys.argv = sys.argv[:1]
    unittest.main()
