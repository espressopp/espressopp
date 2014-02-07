"""
   This script is an example of pressure tensor calculation layerwise according to the
 Irvin Kirwood method.
 Initial configuration file is 'lennard_jones.xyz' (equilibrated lennard-jones fluid).
"""

import espresso
import mpi4py.MPI as MPI

# skin for Verlet list
skin = 0.3
# LJ cutoff
rc   = 2.5
# integration step
dt   = 0.005

# read a configuration from a file
pid, type, xpos, ypos, zpos, xvel, yvel, zvel, Lx, Ly, Lz = \
        espresso.tools.readxyz('lennard_jones.xyz')
# number of particles
NPart              = len(xpos)
# system box size
box                = (Lx, Ly, Lz)
# create a basic system
system             = espresso.System()
# specify a random number generator
system.rng         = espresso.esutil.RNG()
# use orthorhombic periodic boundary conditions
system.bc          = espresso.bc.OrthorhombicBC(system.rng, box)
system.skin        = skin

comm = MPI.COMM_WORLD

nodeGrid = espresso.tools.decomp.nodeGrid(comm.size)
# calculate a 3D subgrid to speed up verlet list builds and communication
cellGrid           = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
# create a domain decomposition particle storage with the specified nodeGrid and cellGrid
system.storage     = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

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
  part = [pid[i], type[i], espresso.Real3D(xpos[i], ypos[i], zpos[i]), espresso.Real3D(xvel[i], yvel[i], zvel[i])]
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
vl      = espresso.VerletList(system, cutoff = rc)
potLJ   = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc)
interLJ = espresso.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# use a velocity Verlet integration scheme
integrator     = espresso.integrator.VelocityVerlet(system)
# set the integration step
integrator.dt  = dt

integrator.run(1)

print "system setup finished"

print 'Calculating pressure...'

n = int(10)         # 10 layers in z direction
h0 = Lz / float(n)  # z coordinate of initial layer
dh = 3.             # area around the layer for kinetic part of the pressure tensor

pressure_tensor    = espresso.analysis.PressureTensor(system)
pressure_tensor_l  = espresso.analysis.PressureTensorLayer(system, h0, dh)
pressure_tensor_ml = espresso.analysis.PressureTensorMultiLayer(system, n, dh)

n_measurements = 10 # result will be averaged over 10 measurements

print 'result will be averaged over ', n_measurements, ' measurements'

pij_layers1 = []
pij_layers2 = []
Pijtot = espresso.Tensor(0.0)

for i in range(n_measurements):
  integrator.run(10)
  print 'measurement Nr: %d of %d' % ( i+1, n_measurements )
  
  # compute the total pressure tensor
  Pijtot += espresso.Tensor( pressure_tensor.compute() )
  
  # layerwise
  pij_aux = pressure_tensor_ml.compute()
  
  for j in range(n):
    pressure_tensor_l.h0 = j * h0
    if(j>= len( pij_layers1 ) ):
      pij_layers1.append( espresso.Tensor( pressure_tensor_l.compute() ) )
      pij_layers2.append( pij_aux[j] )
    else:
      pij_layers1[j] += espresso.Tensor( pressure_tensor_l.compute() )
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
