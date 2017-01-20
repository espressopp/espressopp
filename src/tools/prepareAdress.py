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


"""

***************************************
prepareAdress - setup AdResS simulation
***************************************

Auxiliary python functions for preparation of an Adress Simulation based on a configuration from an all-atomistic simulation.

If one uses a configuration file from an all-atomistic simulation as start configuration for an AdResS simulation,
the particles are probably all located inside the simulation box. However, in AdResS only the coarse-grained center-of-mass
particles have to be in the box, the atomistic particles of the coarse grained might be outside around their CoM CG particle. When in the start configuration
atomistic particles belonging to a molecule are folded such that some of the atoms are on the one side of the box while the others are folded
to the other side the calculation of the center of mass goes wrong and the simulation will be incorrect. This script ensures a proper
center of mass calculation and a proper folding and configuration for the AdResS simulation by simply putting the CG particle in one of the atoms (AdressSetCG) first.
Then the molecules will be put together properly afterwards when calling AdressDecomp.

"""

import espressopp

def AdressSetCG(apm, pidCG, allParticlesAT):
    cmp = [0,0,0]
    pos = (allParticlesAT[pidCG*apm])[1]
    for i in xrange(3):
	cmp[i] += pos[i]
    return cmp

def AdressDecomp(system, integrator):
    system.storage.decompose()
    integrator.run(0)
    system.storage.decompose()
    integrator.run(0)
