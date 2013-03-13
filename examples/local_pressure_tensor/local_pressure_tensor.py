"""
   This script is an example of pressure tensor calculation layerwise according to the
 Irvin Kirwood method.
 Layers should be
 perpendicular to the z-direction of the simulation box. Information about simulation
 box and particles will be read from configuration file 'lennard_jones.xyz'. It is 
 already equilibrated lennard-jones fluid.
"""

import espresso
import MPI

# define some global values
skin = 0.3
# this is the cutoff of our LJ-interaction
rc   = 2.5
# integration step
dt   = 0.005

# read the configuration from a file
pid, type, xpos, ypos, zpos, xvel, yvel, zvel, Lx, Ly, Lz = \
        espresso.tools.readxyz('lennard_jones.xyz')
# we can get the number of particles of the system from the length of the pid-list
NPart              = len(xpos)
# get the box size from the file
box                = (Lx, Ly, Lz)
# create the basic system
system             = espresso.System()
# we always have to specify a random number generator (even if we do not need it in this example)
system.rng         = espresso.esutil.RNG()
# use orthorhombic periodic boundary conditions
system.bc          = espresso.bc.OrthorhombicBC(system.rng, box)
# we also have to provide a skin for things like Verlet-Lists
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
  #part = [pid[i], type[i], espresso.Real3D(xpos[i], ypos[i], zpos[i]+Lz/2.), espresso.Real3D(xvel[i], yvel[i], zvel[i])]
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

# setup the analysis for the pressure tensor
pressure_tensor = espresso.analysis.PressureTensor(system)

print 'Calculating pressure...'

n = int(10)         # we will calculate pressure in 10 layers in z direction
z0 = Lz / float(n)  # z coordinate of initial layer
dz = 3.             # area around the layer where the pressure will be calculated
n_measurements = 10 # result will be averaged over 100 maesurements

print 'result will be averaged over ', n_measurements, ' measurements'

pij_layers1 = []
pij_layers2 = []
Pijtot = espresso.Tensor(0.0)
for i in range(n_measurements):
  integrator.run(10)
  print 'measurement Nr:', (i+1), 'of', n_measurements
  
  # compute the tensor for whole box
  Pijtot += pressure_tensor.compute()
  
  # layerwise
  pij_aux = pressure_tensor.compute(n, dz)
  
  for j in range(n):
    if(i==0):
      pij_layers1.append( pressure_tensor.compute( j * z0, dz) )
      pij_layers2.append( pij_aux[j] )
    else:
      pij_layers1[j] += pressure_tensor.compute( j * z0, dz)
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

print '\npressure layerwise; this is the function where one should set a Z coordinate of the plane (float number) \n \
where pressure tensor will be calculated. (Advantage: one can set the position precisely)'
      
print 'layer number     z coord of layer        pressure tensor'
for i in range(n):
  print ('%4d           %7.3f              ' % (i, i * z0)) , pij_layers1[i]

print '\npressure layerwise; this is the function where one should set a number of layers (integer number N). \n \
Lz will be devided by N and then pressure tensor will be calculated in each layer. (Advantage: it is faster)'

print 'layer number     z coord of layer        pressure tensor'
for i in range(n):
  print ('%4d           %7.3f              ' % (i, i * z0)) , pij_layers2[i]
  
print 'done'
print 'both functions should give the same result'
