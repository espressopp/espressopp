#!/usr/bin/env python2
#  Copyright (C) 2016-2017(H)
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
# DEMONSTRATION OF THE LATTICE-BOLTZMANN SIMULATION			  #	
###########################################################################

import espressopp
from espressopp import Int3D
from espressopp import Real3D

L = 20
tempT = 1.
num_particles = 100
# create default Lennard Jones (WCA) system with 0 particles and cubic box
system, integrator = espressopp.standard_system.LennardJones(num_particles, box=(L, L, L), temperature=tempT)
integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = 0.005

# define a LB grid
nodeGrid=espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
integrator.addExtension(lb)

# set parameters of LB fluid (LJ units)
lb.lbTemp = tempT       # desired temperature
lb.nSteps = 10          # time step contrast between LB and MD (t_{LB} / t_{MD})
lb.visc_b = 3.          # bulk viscosity of LB fluid
lb.visc_s = 3.          # shear viscosity of LB fluid
lb.profStep = 5000      # time profiling frequency

# initialize populations
initDen = 1.
initVel = Real3D (0.)
#initPop = espressopp.integrator.LBInitPopWave(system,lb)   # sin-like in z-dir
initPop = espressopp.integrator.LBInitPopUniform(system,lb) # uniform
initPop.createDenVel( initDen, initVel )

# screen output
lboutputScreen = espressopp.analysis.LBOutputScreen(system,lb)
outScreen=espressopp.integrator.ExtAnalyze(lboutputScreen,lb.profStep)
integrator.addExtension(outScreen)

# run simulations
for k in range (3):
    integrator.run(50000)
    s = str(integrator.step)
    # output md configuration
    mdoutput = 'dump.' + s + '.xyz'
    espressopp.tools.fastwritexyz(mdoutput, system)
    # output LB configuration
    lb.keepLBDump()     # flag to keep previously saved LB state
    lb.saveLBConf()     # saves current state of the LB fluid
