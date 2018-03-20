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

system, integrator = espressopp.standard_system.LennardJones(0,(10,10,10))
system.storage.addParticles([[1,0,espressopp.Real3D(5,5,5)]],'id','type','pos')
system.storage.addParticles([[2,0,espressopp.Real3D(5,6,5)]],'id','type','pos')
system.storage.addParticles([[3,1,espressopp.Real3D(6,5,5)]],'id','type','pos')
system.storage.addParticles([[4,1,espressopp.Real3D(6,6,5)]],'id','type','pos')
system.storage.decompose()

Epot = espressopp.analysis.EnergyPot(system)
print 'Epot = ',Epot.compute()

intLJ = system.getInteraction(0)
pot00 = intLJ.getPotential(0,0)
intLJ.setPotential(1,1,pot00)
intLJ.setPotential(1,0,pot00)
intLJ.setPotential(0,1,pot00)

print 'Epot = ', Epot.compute()

