#!/usr/bin/env python2
#  Copyright (C) 2015-2017(H)
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

#################################################################################
#                                                                               #
#  ESPResSo++ Python script for a L-J standadrd system featuring force-capping  # 
#                                                                               #
#################################################################################


import espressopp

# try to run this script with force_capping switched on and off
forceCapping = True

# specify number of particles
num_particles = 200
rho = 0.85
L   = pow(num_particles/rho, 1.0/3.0)

# setup random LennardJones system
# this system is likely to explode on integration, because particles can strongly overlap
system, integrator = espressopp.standard_system.LennardJones(num_particles, box=(L, L, L), temperature=1.0)

# choose a smaller timestep
integrator.dt = 0.0001

if forceCapping:
  max_force = 100000.0
  # define force capping extension
  capForce = espressopp.integrator.CapForce(system, max_force)
  # and add it to the integrator
  integrator.addExtension(capForce)

espressopp.tools.analyse.info(system, integrator)

sock = espressopp.tools.vmd.connect(system)
for i in range(1000):
  # make 10 Velocity-Verlet integration steps
  integrator.run(10)
  # print system information
  espressopp.tools.analyse.info(system, integrator)
  # update postions in VMD
  espressopp.tools.vmd.imd_positions(system, sock)
  #switch off force capping after some time
  if forceCapping and i>100:
    capForce.disconnect()
    forceCapping = False
    print "switching off force capping"

