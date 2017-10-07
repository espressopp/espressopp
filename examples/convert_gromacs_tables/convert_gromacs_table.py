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
#  ESPResSo++                                                             #
#  Test script for Converting GROMACS tabulated file                      #
#                                                                         #
###########################################################################


import sys
import time
import os
import espressopp
import mpi4py.MPI as MPI
import math
import logging
import os
from espressopp import Real3D, Int3D
from espressopp.tools import gromacs


# Input values for system
N = 10                                    # box size
size  = (float(N), float(N), float(N))
numParticles = N**3                       # number of particles
nsteps  = 1000                            # number of steps
cutoff  = 2.5                             # cutoff for LJ potential
tabfile = "pot-lj-esp.tab"                # filename for tabulated potential
skin    = 0.3                             # skin for Verlet lists
spline  = 2                               # interpolation spline type

# parameters to convert GROMACS tabulated potential file
filein = "table6-12.xvg"                  # gromacs tabulated file to be converted
fileout = "pot-lj-gro.tab"                # filename of output file
sigma = 1.0
epsilon = 1.0
c6 = 4.0
c12 = 4.0

files = [tabfile, fileout]                # run simulation on these files


######################################################################
##  IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE  ##
######################################################################

print '\n-- GROMACS Tabulated File Conversion Test -- \n'
print 'Steps: %3s' % nsteps
print 'Particles: %3s' % numParticles
print 'Cutoff: %3s' % cutoff


# writes the tabulated potential file
def writeTabFile(pot, name, N, low=0.0, high=2.5, body=2):
    outfile = open(name, "w")
    delta = (high - low) / (N - 1)
     
    for i in range(N):
        r = low + i * delta
        energy = pot.computeEnergy(r)
        if body == 2:# this is for 2-body potentials
            force = pot.computeForce(Real3D(r, 0.0, 0.0))[0]
        else: # this is for 3- and 4-body potentials
            force = pot.computeForce(r)
        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))
     
    outfile.close()

# write the espressopp++ tabulated file for a LJ potential
print 'Generating potential file ... (%2s)' % tabfile
potLJ = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift=0.0, cutoff=cutoff)
writeTabFile(potLJ, tabfile, N=1500, low=0.01, high=potLJ.cutoff)


# convert gromacs tabulated file to espressopp++ format
print 'Converting GROMACS file to ESPResSo++ file ... (%2s -> %2s)' % (filein, fileout)
gromacs.convertTable(filein, fileout, sigma, epsilon, c6, c12)

#exit() # exit if you just want to convert a file

# compute the number of cells on each node
def calcNumberCells(size, nodes, cutoff):
    ncells = 1
    while size / (ncells * nodes) >= cutoff:
       ncells = ncells + 1
    return ncells - 1


#start_time = time.clock()

# run simulation for all tabulated potential files
for potfile in files:
    print '\nUsing file: %0s'% potfile
        
    # set up system
    system = espressopp.System()
    system.rng  = espressopp.esutil.RNG()
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
    system.skin = skin
        
    comm = MPI.COMM_WORLD
        
    nodeGrid = Int3D(1, 1, comm.size)
    cellGrid = Int3D(
        calcNumberCells(size[0], nodeGrid[0], cutoff),
        calcNumberCells(size[1], nodeGrid[1], cutoff),
        calcNumberCells(size[2], nodeGrid[2], cutoff)
        )
    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
        
    pid = 0
        
    for i in range(N):
        for j in range(N):
            for k in range(N):
                m = (i + 2*j + 3*k) % 11
                r = 0.45 + m * 0.01
                x = (i + r) / N * size[0]
                y = (j + r) / N * size[1]
                z = (k + r) / N * size[2]
                    
                x = 1.0 * i
                y = 1.0 * j
                z = 1.0 * k
                    
                system.storage.addParticle(pid, Real3D(x, y, z))
                    
                # not yet: dd.setVelocity(id, (1.0, 0.0, 0.0))
                pid = pid + 1
        
        
    system.storage.decompose()
        
    # integrator
    integrator = espressopp.integrator.VelocityVerlet(system)
    integrator.dt = 0.005
        
    # now build Verlet List
    # ATTENTION: you must not add the skin explicitly here
    logging.getLogger("Interpolation").setLevel(logging.INFO)
    vl = espressopp.VerletList(system, cutoff = cutoff + system.skin)
        
    potTab = espressopp.interaction.Tabulated(itype=spline, filename=potfile, cutoff=cutoff)
    # ATTENTION: auto shift was enabled
    interTab = espressopp.interaction.VerletListTabulated(vl)
    interTab.setPotential(type1=0, type2=0, potential=potTab)
    system.addInteraction(interTab)
        
        
    temp = espressopp.analysis.Temperature(system)
    press = espressopp.analysis.Pressure(system)
        
    temperature = temp.compute()
    p = press.compute()
    Ek = 0.5 * temperature * (3 * numParticles)
    Ep = interTab.computeEnergy()
        
        
    print 'Start %5s: tot energy = %10.3f pot = %10.3f kin = %10.3f temp = %10.3f p = %10.3f' \
        % ("", Ek + Ep, Ep, Ek, temperature, p)
        
    # langevin thermostat
    langevin = espressopp.integrator.LangevinThermostat(system)
    integrator.addExtension(langevin)
    langevin.gamma = 1.0
    langevin.temperature = 1.0
        
    integrator.run(nsteps)
    temperature = temp.compute()
    p = press.compute()
    Ek = 0.5 * temperature * (3 * numParticles)
    Ep = interTab.computeEnergy()
    print 'Step %6d: tot energy = %10.3f pot = %10.3f kin = %10.3f temp = %10.3f p = %10.3f' % \
            (nsteps, Ek + Ep, Ep, Ek, temperature, p)
       
    os.system('rm '+potfile) # remove file

print '\nDone.'
