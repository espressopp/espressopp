"""
We just want to read in a particle configuration of an equilibrated
lennard-jones fluid and calculate the pressure tensor in slabs along
the z-direction of the simulation box. This can be done by manually
specifying a slab node grid for the domain decomposition storage and
calculating the pressure tensor for each node (cpu) locally.
"""

import espresso

# define some global values
skin = 0.3
# this is the cutoff of our LJ-interaction
rc   = 2.5
# integration step
dt   = 0.005

# read the configuration from a file
pid, type, xpos, ypos, zpos, xvel, yvel, zvel, Lx, Ly, Lz = espresso.tools.readxyz('lennard_jones.xyz')
# we can get the number of particles of the system from the length of the pid-list
NPart              = len(pid)
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
# get the number of CPUs to use
NCPUs              = espresso.MPI.COMM_WORLD.size
# manually specify the node grid e.g. slabs in z-direction
# (we want to calculate the pressure tensor in slabs)
nodeGrid           = (1,1,NCPUs)
# calculate a 3D subgrid to speed up verlet list builds and communication
cellGrid           = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
# create a domain decomposition particle storage with the specified nodeGrid and cellGrid
system.storage     = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

print "NCPUs              = ", NCPUs
print "nodeGrid           = ", nodeGrid
print "cellGrid           = ", cellGrid
print "NPart              = ", NPart
print "skin               = ", skin
print "rc                 = ", rc
print "dt                 = ", dt

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
potLJ   = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=0.0)
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
Pxy = espresso.analysis.PressureTensor(system)
# compute the tensor locally and return each nodes values in a python list
LocalPxy = Pxy.computeLocal()

for i in range(len(LocalPxy)):
  nodePxy = LocalPxy[i]
  # the volume box of the nodes domain is also returned by the command
  volumeBox   = nodePxy[0]
  # t is the local tensor
  t    = nodePxy[1]
  xmin = volumeBox[0]
  xmax = volumeBox[1]
  ymin = volumeBox[2]
  ymax = volumeBox[3]
  zmin = volumeBox[4]
  zmax = volumeBox[5]
  print "volume(%f %f %f %f %f %f): Pxx= %f Pyy= %f Pzz= %f Pxy= %f Pxz= %f Pyz= %f" % (xmin, xmax, ymin, ymax, zmin, zmax, t[0], t[1], t[2], t[3], t[4], t[5])
