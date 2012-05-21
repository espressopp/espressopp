'''
  Here we check storage.removeParticle(id)
'''

import sys
import time
import espresso
import MPI
from espresso.tools.convert import lammps
from espresso.tools import decomp

'''
if with_removing is false the 27 particles from file 'particles27.start' will be added and standard
system with Lennard Jones potential will be run.
if with_removing is true then 4 random particles will be added, then 27 particles from file
'particles27.start' and then more 6 random particles will be added. In this case before running
the simulation all the random particles will be removed. Thereby at the begin of simulation we will 
have 27 particles from 'particles27.start' only.
Thus one can compare two identical systems (They should be identical).

Resulting file 'result.dat' consist of 6 column:
timestep, temperature, pressure, total energy, potential energy, kinetic energy
'''
with_removing = False

# cutoff, skin, integration timestep
rc = pow(2.0, 1./6.)
skin = 0.3
timestep = 0.0045

# reading the lammps file with 27 particles
file = 'particles27.start'
# coordinates, size of the box, velosities
x, y, z, Lx, Ly, Lz, vx, vy, vz = lammps.read(file)

num_particles = len(x)

####################################
### SETTING THE USUAL PARAMETERS ###
####################################
sys.stdout.write('Setting up simulation ...\n')
density = num_particles / (Lx * Ly * Lz)
box = (Lx, Ly, Lz)
system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, box)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(box, nodeGrid, rc, skin)

system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
####################################

# add particles to the system and then decompose
props = ['id', 'type', 'pos', 'v', 'mass']
new_particles = []

if with_removing:
  
  for i in range (num_particles+10):
    pid = i+1
    type = 0
    mass = 1.0
    if(i<4 or i>=31):
      # adding 4 random particles at the begin and 6 random particles at the end
      pos  = system.bc.getRandomPos()
      vel  = espresso.Real3D(0.0, 0.0, 0.0)
    else:
      # adding fixed particles from file
      pos  = espresso.Real3D(x[i-4],y[i-4],z[i-4])
      vel  = espresso.Real3D(vx[i-4],vy[i-4],vz[i-4])
    part = [pid, type, pos, vel, mass]
    new_particles.append(part)
  
else:
  
  # adding just 27 particles from the file
  for i in range (num_particles):
    pid = i + 1
    type = 0
    mass = 1.0
    pos  = espresso.Real3D(x[i],y[i],z[i])
    vel  = espresso.Real3D(vx[i],vy[i],vz[i])
    part = [pid, type, pos, vel, mass]
    new_particles.append(part)
    
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

# all particles interact via a LJ interaction (use Verlet lists)
vl = espresso.VerletList(system, cutoff=rc+system.skin)
potLJ = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, cutoff=rc, shift=True)
interLJ = espresso.interaction.VerletListLennardJones(vl)
interLJ.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interLJ)

# setup integrator
integrator = espresso.integrator.VelocityVerlet(system)
integrator.dt = timestep

print ''
print 'number of particles =', num_particles
print 'box: ', box
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)

formatPrint = '%5d %8.4f %10.5f %12.3f %12.3f %12.3f'

T = temperature.compute()
P = pressure.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interLJ.computeEnergy()
print(' step     T          P       etotal     epotential    ekinetic\n')
print(formatPrint % (0, T, P, Ek + Ep, Ep, Ek))

if with_removing:
  for i in range (4):
    system.storage.removeParticle(i+1)
  for i in range (6):
    system.storage.removeParticle(i+1+31)

new_num_of_p = espresso.analysis.NPart(system).compute()
print 'new number of particles: ', new_num_of_p
  
print '\n Equlibration!'
# initial equilibration
for i in range (10):
  integrator.run(10000)
  step = 10000 * (i+1)
  T = temperature.compute()
  P = pressure.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interLJ.computeEnergy()
  print(formatPrint % (step, T, P, Ek + Ep, Ep, Ek))

# here we start the measurement
print '\n Measurements!'

# resulting file
fmt2 = ' %10d %6.3f %6.3f %15.3f %15.3f %15.3f\n'
if with_removing:
  nameFileK = 'resultWith.dat'
else:
  nameFileK = 'resultWithout.dat'
resFileK = open (nameFileK, 'a')

for i in range (20):
  
  for j in range (100):
    integrator.run(100)
    instep = 100*j+i*10000
    T = temperature.compute()
    P = pressure.compute()
    Ek = 0.5 * T * (3 * num_particles)
    Ep = interLJ.computeEnergy()
    resFileK.write(fmt2 % ( instep, T, P, (Ek+Ep), Ep, Ek))
  
  step = (i+1)*10000
  T = temperature.compute()
  P = pressure.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interLJ.computeEnergy()
  
  print(formatPrint % (step, T, P, Ek + Ep, Ep, Ek))
  
resFileK.close()