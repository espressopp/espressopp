#!/usr/bin/env python2
#
#  Copyright (C) 2013-2017(H)
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
# 
# -*- coding: utf-8 -*-

import espressopp
from espressopp import Real3D

d = 0.85
Nchains = 10

Mmonomers = 10
N = Nchains * Mmonomers
L = pow(N/d, 1.0/3)

system, integrator = espressopp.standard_system.PolymerMelt(Nchains, Mmonomers,(10,10,10), dt = 0.005, temperature=1.0)


print "starting warmup"
org_dt = integrator.dt
pot = system.getInteraction(0).getPotential(0,0)
print pot
print "Nint = ", system.getNumberOfInteractions()
final_sigma = pot.sigma
final_epsilon = pot.epsilon
print "sigma=",pot.sigma, "epsilon=",pot.epsilon
maxParticleID = int(espressopp.analysis.MaxPID(system).compute())
N = 1
number = 50

for k in range(number):
  if k < 10:
    continue
 
  force_capping = espressopp.integrator.CapForce(system, 1000000.0/number*k)
  integrator.addExtension(force_capping)

  pot.sigma = final_sigma/number*k
  pot.epsilon = final_epsilon/number*k

  integrator.dt = 0.0001
  espressopp.tools.analyse.info(system, integrator)
  integrator.run(N)
  espressopp.tools.analyse.info(system, integrator)

integrator.dt = org_dt
pot.sigma = final_sigma
pot.epsilon = final_epsilon
force_capping.disconnect()
  
for k in range(10):
  integrator.run(70)
  espressopp.tools.analyse.info(system, integrator)
integrator.step = 0

print "warmup finished"

for k in range(10):
  integrator.run(100)
  espressopp.tools.analyse.info(system, integrator)

