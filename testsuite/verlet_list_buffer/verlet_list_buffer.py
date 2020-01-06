#!/usr/bin/env python2

import unittest
import espressopp
from espressopp.tools import readxyz
import time

def generate_vl(useBuffers):
  print 'VERLET LIST {}USING BUFFERS'.format('NOT ' if not useBuffers else '')
  nsteps      = 1
  isteps      = 1000
  rc          = 2.5
  skin        = 0.4
  timestep    = 0.005
  dt          = 0.005
  epsilon     = 1.0
  sigma       = 1.0

  # set temperature to None for NVE-simulations
  temperature = 1.0

  xyz_file = "lennard_jones_fluid_10000.xyz"
  pid, type, xpos, ypos, zpos, xvel, yvel, zvel, Lx, Ly, Lz = readxyz(xyz_file)
  box = (Lx, Ly, Lz)
  num_particles = len(pid)

  system, integrator = espressopp.standard_system.Default(box=box, rc=rc, skin=skin, dt=timestep, temperature=temperature)

  props = ['id', 'type', 'mass', 'pos', 'v']
  new_particles = []
  for i in range(num_particles):
    part = [i + 1, 0, 1.0, espressopp.Real3D(xpos[i], ypos[i], zpos[i]), espressopp.Real3D(xvel[i], yvel[i], zvel[i])]
    new_particles.append(part)
    if i % 1000 == 0:
      system.storage.addParticles(new_particles, *props)
      system.storage.decompose()
      new_particles = []
  system.storage.addParticles(new_particles, *props)
  system.storage.decompose()

  # Lennard-Jones with Verlet list
  vl      = espressopp.VerletList(system, cutoff = rc, useBuffers = useBuffers)
  potLJ   = espressopp.interaction.LennardJones(epsilon=epsilon, sigma=sigma, cutoff=rc, shift=0)
  interLJ = espressopp.interaction.VerletListLennardJones(vl)
  interLJ.setPotential(type1=0, type2=0, potential=potLJ)
  system.addInteraction(interLJ)

  # espressopp.tools.analyse.info(system, integrator)

  espressopp.tools.analyse.info(system, integrator)
  start_time = time.clock()
  for k in range(nsteps):
    integrator.run(isteps)
    espressopp.tools.analyse.info(system, integrator)
  end_time = time.clock()
  espressopp.tools.analyse.final_info(system, integrator, vl, start_time, end_time)

  pairs = sum(vl.getAllPairs(),[])
  return pairs

def sort_pairs(pairs):
  # sort each tuple
  pairs = [(p[1],p[0]) if p[1]<p[0] else p for p in pairs]
  # sort all tuples in list
  pairs = sorted(pairs)
  return pairs

class TestVerletListBuffer(unittest.TestCase):
  def test1vl(self):
    print '-'*70
    pairs1 = sort_pairs(generate_vl(False))
    print '-'*70
    pairs2 = sort_pairs(generate_vl(True))

    # ensure the same pairs are generated
    self.assertEqual(len(pairs1), len(pairs2))
    for i in range(len(pairs1)):
      self.assertEqual(pairs1[i],pairs2[i])

if __name__ == "__main__":
  unittest.main()
