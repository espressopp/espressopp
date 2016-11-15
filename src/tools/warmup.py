#  Copyright (C) 2012,2013,2016
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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


import espressopp
from espressopp import Real3D

def warmup(system, integrator, number = 80):
  """
  Warm up for a system with a density of 0.85.

  The method needs the following parameters:

  * system, integrator
    ESPResSo system which schoul be warmed up and the correspondig integrator e.g.:
    
    >>> system, integrator = espressopp.standard_system.LennardJones(100,(10,10,10))
    
  * number
    number of steps of the warm up 
    
    for a system with a density of 0.85, if it explodes try a higher number
    
  """
  print "starting warmup"

  org_dt = integrator.dt
  pot = system.getInteraction(0).clonePotential(0,0)
  final_sigma = pot.sigma
  final_epsilon = pot.epsilon
  N = 50

  integrator.dt = 0.0001
  force_capping = espressopp.integrator.CapForce(system, 0.0)
  integrator.addExtension(force_capping)

  for k in xrange(11,number):
    force_capping.setAbsCapForce(1000000.0/number*k)
    pot.sigma = final_sigma/number*k
    pot.epsilon = final_epsilon/number*k
    system.getInteraction(0).setPotential(0,0,pot)
    espressopp.tools.analyse.info(system, integrator)
    integrator.run(N)

  integrator.dt = org_dt
  pot.sigma = final_sigma
  pot.epsilon = final_epsilon
  force_capping.disconnect()
  
  for k in xrange(11):
    integrator.run(70)
    espressopp.tools.analyse.info(system, integrator)

  integrator.step = 0

  print "warmup finished"
