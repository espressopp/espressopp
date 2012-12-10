"""
We just want to read in a particle configuration of an equilibrated
lennard-jones fluid and calculate the pressure tensor in slabs along
the z-direction of the simulation box.
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

comm = MPI.COMM_WORLD

nodeGrid = espresso.tools.decomp.nodeGrid(comm.size)
# calculate a 3D subgrid to speed up verlet list builds and communication
cellGrid           = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
# create a domain decomposition particle storage with the specified nodeGrid and cellGrid
system.storage     = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

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
# compute the tensor for whole box
Pijtot = Pxy.compute()
print 'total tensor'
print '   Pxx      Pyy      Pzz      Pxy      Pxz      Pyz'
fmt1 = '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f'
print(fmt1 % (Pijtot[0], Pijtot[1], Pijtot[2], Pijtot[3], Pijtot[4], Pijtot[5]))

# compute the tensor locally and return each nodes values in a python list
n = 30; # we will divide the box in 30 slabs (z direction)
dLz = Lz / float(n)
Pijloc = []
for i in range(n):
  # The rutine .compute(xmin, xmax, ymin, ymax, zmin, zmax) will return the 
  # pressure tensor for defined volume
  Pij = Pxy.compute(0.0, Lx, 0.0, Ly, dLz*i, dLz*(i+1))
  Pijloc.append(Pij)

fmt2 = '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f'
print 'local tensor'
print '  zcoord      Pxx      Pyy      Pzz      Pxy      Pxz      Pyz'
for i in range(n):
  print(fmt2 % ( dLz*(i+0.5), Pijloc[i][0], Pijloc[i][1], Pijloc[i][2], Pijloc[i][3], Pijloc[i][4], Pijloc[i][5]))
