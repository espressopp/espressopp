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


'''
#  This script is an example of calculation of long range interactions (Coulomb interaction) using 
#  the Ewald summation method.
#  
#  The initial data file is 'ini_struct_deserno.dat'. It contains initial information about
#  the system: the box size, particle id, position and charge.
#
#  File 'deserno_ewald.dat' contains results which were obtained by Markus Deserno for exactly
#  the same system. Ewald parameter (alpha = 1.112583061), Cutoff in R space (rspacecutoff = 4.9)
#  cutoff in K space (kspacecutoff = 30). It compares the results of Markus Deserno with Espresso++
#  implementation. The forces which affect the particles, total energy and the differences between
#  Deserno's results and Espresso++ implementation are printed at the end of the script.
'''

# this is an auxiliary function. It reads the results of Deserno from "deserno_ewald.dat"
def readingDesernoForcesFile():
  # forces x,y,z
  fx, fy, fz = [], [], []
  # energy
  energy = 0.0
  
  # reading the general information
  file = open("deserno_ewald.dat")
  i = 0
  for line in file:
    # energy
    if i==6:
      tmp = line.split()
      energy = float(tmp[0])
    
    # forces
    if i>=9:
      line = line.replace('{','').replace('}','')
      tmp = line.split()
      fx.append(float(tmp[0]))
      fy.append(float(tmp[1]))
      fz.append(float(tmp[2]))
      
    i=i+1
  
  return energy, fx, fy, fz
# end of the function readingDesernoForcesFile


# The script itself
import sys
import mpi4py.MPI as MPI
import espressopp

from espressopp import Real3D
from espressopp.tools import espresso_old

# reading the particle coordinates, charges and box size from old espressopp data file
# file 'ini_struct_deserno.dat' contains the data we need
print "Reading system data:"
Lx, Ly, Lz, x, y, z, type, q, vx,vy,vz,fx,fy,fz,bondpairs = espresso_old.read('ini_struct_deserno.dat')

# creating the system box
box = (Lx, Ly, Lz)
print "System box size:", box
# number of particles
num_particles = len(x)
print "Number of particles = ", num_particles
print "The first particle has coordinates", x[0], y[0], z[0]

'''
#  Ewald method suppose to calculate electrostatic interaction dividing it into R space and
#  K space part
#  
#  alpha - Ewald parameter
#  rspacecutoff - the cutoff in real space
#  kspacecutoff - the cutoff in reciprocal space   
'''
alpha          = 1.112583061
rspacecutoff   = 4.9
kspacecutoff   = 30

# seting the skin for Verlet list (it is not important here)
skin           = 0.09

# Coulomb prefactor parameters
bjerrumlength  = 1.0
temperature    = 1.0
coulomb_prefactor = bjerrumlength * temperature

nodeGrid       = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size,box,rspacecutoff,skin)
cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rspacecutoff, skin)
system         = espressopp.System()
system.rng     = espressopp.esutil.RNG()
system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# Adding the particles
props = ['id', 'pos', 'type', 'q']
new_particles = []
for i in range(0, num_particles):
  part = [ i, Real3D(x[i], y[i], z[i]), type[i], q[i] ]
  new_particles.append(part)
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

## potential and interaction ##

# setting the Verlet list
vl = espressopp.VerletList(system, rspacecutoff+skin)

# the R space part of electrostatic interaction according to the Ewald method
'''
  Creating the Coulomb potential which is responsible for the R space part according to the
  Ewald method.
  It is based on the Coulomb prefactor (coulomb_prefactor), Ewald parameter (alpha),
  and the cutoff in R space (rspacecutoff)
'''
coulombR_pot = espressopp.interaction.CoulombRSpace(coulomb_prefactor, alpha, rspacecutoff)
# creating the interaction based on the Verlet list
coulombR_int = espressopp.interaction.VerletListCoulombRSpace(vl)
# setting the potential for the interaction between particles of type 0 and 0
coulombR_int.setPotential(type1=0, type2=0, potential = coulombR_pot)
# adding the interaction to the system
system.addInteraction(coulombR_int)

# k space part of electrostatic interaction according to the Ewald method
'''
  Creating the Coulomb potential which is responsible for the K space part according to the
  Ewald method.
  It is based on the system information (system), Coulomb prefactor (coulomb_prefactor),
  Ewald parameter (alpha), and the cutoff in K space (kspacecutoff)
'''
ewaldK_pot = espressopp.interaction.CoulombKSpaceEwald(system, coulomb_prefactor, alpha, kspacecutoff)
# creating the interaction based on the Cell list for all particle interaction and potential in K space
ewaldK_int = espressopp.interaction.CellListCoulombKSpaceEwald(system.storage, ewaldK_pot)
# adding the interaction to the system
system.addInteraction(ewaldK_int)

# creating the integrator which based on the Verlet algorithm
integrator    = espressopp.integrator.VelocityVerlet(system)
# seting the time step (it is not important here)
integrator.dt = 0.0001

# nothing will be changed in system, just forces will be calculated ones
integrator.run(0)

# reading Deserno results (energy and forces)
energy_Deserno, forceX_Deserno, forceY_Deserno, forceZ_Deserno = readingDesernoForcesFile()

# printing the particle id, force (x,y,z), and force difference (x,y,z)
format0 = '\n %45s %105s \n'
print (format0 % ('forces', 'the difference between Deserno\'s result and forces by Espresso++'))
format1 = '%3s %20s %20s %20s %10s %20s %25s %25s\n'
print (format1 % ('id', 'fx', 'fy', 'fz', ' ', 'dfx', 'dfy', 'dfz'))
format2 = '%3d %3s %3.17f %3s %3.17f %3s %3.17f %10s %3.17f %3s %3.17f %3s %3.17f'
for j in range(0, num_particles):
  print (format2 % (j, ' ', \
                       system.storage.getParticle(j).f.x, ' ', \
                       system.storage.getParticle(j).f.y, ' ', \
                       system.storage.getParticle(j).f.z, \
                       ' ', \
                       abs(system.storage.getParticle(j).f.x-forceX_Deserno[j]), ' ', \
                       abs(system.storage.getParticle(j).f.y-forceY_Deserno[j]), ' ', \
                       abs(system.storage.getParticle(j).f.z-forceZ_Deserno[j])) )
  
# calculating the R space part of electrostatic energy
enR = coulombR_int.computeEnergy()
# calculating the K space part of electrostatic energy
enK = ewaldK_int.computeEnergy()
# total energy
enTot = enR + enK

# printing the total energy and the difference with Deserno results
print '\nTotal energy: %5.16f;         The difference in energy (Deserno\'s result, Espresso++): %5.16f\n' % (enTot, enTot-energy_Deserno)

sys.exit()
