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

###########################################################################			
#                                                                         #
#  ESPResSo++ Python script for an example of pressure tensor calculation #
#  layerwise according to the  Irvin Kirwood method			  #
#                                                                         #
###########################################################################

"""
 Initial configuration file is 'lennard_jones.xyz' (equilibrated lennard-jones fluid).
"""

import espressopp
import mpi4py.MPI as MPI

# skin for Verlet list
skin = 0.3
# LJ cutoff
rc   = 2.5
# integration step
dt   = 0.005

# read a configuration from a file
pid, type, xpos, ypos, zpos, xvel, yvel, zvel, Lx, Ly, Lz = \
        espressopp.tools.readxyz('lennard_jones.xyz')
# number of particles
NPart              = len(xpos)
# system box size
box                = (Lx, Ly, Lz)
# create a basic system
system             = espressopp.System()
# specify a random number generator
system.rng         = espressopp.esutil.RNG()
# use orthorhombic periodic boundary conditions
system.bc          = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin        = skin

comm = MPI.COMM_WORLD

nodeGrid = espressopp.tools.decomp.nodeGrid(comm.size,box,rc,skin)
# calculate a 3D subgrid to speed up verlet list builds and communication
cellGrid           = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
# create a domain decomposition particle storage with the specified nodeGrid and cellGrid
system.storage     = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

print "number of particles = ", NPart
print "box                 = ", box
print "nodeGrid            = ", nodeGrid
print "cellGrid            = ", cellGrid
print "skin                = ", skin
print "cutoff              = ", rc
print "timestep            = ", dt

print "setting up system ..."
# add the particles from the file to the storage of the system
properties = ['id', 'type', 'pos', 'v']
particles  = []
for i in range(NPart):
  part = [pid[i], type[i], espressopp.Real3D(xpos[i], ypos[i], zpos[i]), espressopp.Real3D(xvel[i], yvel[i], zvel[i])]
  particles.append(part)
  # add particles in chunks of 1000 particles, this is faster than adding each single particle
  if i % 1000 == 0:
    system.storage.addParticles(particles, *properties)
    # distribute particles to their nodes
    system.storage.decompose()
    particles = []
system.storage.addParticles(particles, *properties)
system.storage.decompose()

# setup the Lennard-Jones interaction, we use Verlet-List to loop over all interactions
vl      = espressopp.VerletList(system, cutoff = rc)
potLJ   = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc)
interLJ = espressopp.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# use a velocity Verlet integration scheme
integrator     = espressopp.integrator.VelocityVerlet(system)
# set the integration step
integrator.dt  = dt

integrator.run(1)

print "system setup finished"

print 'Calculating pressure...'

n = int(10)         # 10 layers in z direction
h0 = Lz / float(n)  # z coordinate of initial layer
dh = 3.             # area around the layer for kinetic part of the pressure tensor

pressure_tensor    = espressopp.analysis.PressureTensor(system)
pressure_tensor_l  = espressopp.analysis.PressureTensorLayer(system, h0, dh)
pressure_tensor_ml = espressopp.analysis.PressureTensorMultiLayer(system, n, dh)

n_measurements = 10 # result will be averaged over 10 measurements

print 'result will be averaged over ', n_measurements, ' measurements'

pij_layers1 = []
pij_layers2 = []
Pijtot = espressopp.Tensor(0.0)

for i in range(n_measurements):
  integrator.run(10)
  print 'measurement Nr: %d of %d' % ( i+1, n_measurements )
  
  # compute the total pressure tensor
  Pijtot += espressopp.Tensor( pressure_tensor.compute() )
  
  # layerwise
  pij_aux = pressure_tensor_ml.compute()
  
  for j in range(n):
    pressure_tensor_l.h0 = j * h0
    if(j>= len( pij_layers1 ) ):
      pij_layers1.append( espressopp.Tensor( pressure_tensor_l.compute() ) )
      pij_layers2.append( pij_aux[j] )
    else:
      pij_layers1[j] += espressopp.Tensor( pressure_tensor_l.compute() )
      pij_layers2[j] += pij_aux[j]
  

# averaging
Pijtot /= float(n_measurements)
for i in range(n):
  pij_layers1[i] /= float(n_measurements)
  pij_layers2[i] /= float(n_measurements)

print '\ntotal pressure tensor'
print '   Pxx      Pyy      Pzz      Pxy      Pxz      Pyz'
fmt1 = '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f'
print(fmt1 % (Pijtot[0], Pijtot[1], Pijtot[2], Pijtot[3], Pijtot[4], Pijtot[5]))

print '\nPressure tensor by PressureTensorLayer (caculated for each layer separatelly).'
      
print 'layer number     z coord of layer        pressure tensor'
for i in range(n):
  print ('%4d           %7.3f              ' % (i, i * h0)) , pij_layers1[i]

print '\nPressure tensor by PressureTensorMultiLayer (caculated for each layer at once).'

print 'layer number     z coord of layer        pressure tensor'
for i in range(n):
  print ('%4d           %7.3f              ' % (i, i * h0)) , pij_layers2[i]
  
print 'done'
print 'both functions should give the same result'
