'''
#  This script is an example of calculation of long range interactions (Coulomb interaction) using 
#  the Ewald summation method.
'''
import sys
import MPI
import espresso

from espresso import Real3D
from espresso.tools.convert import espresso_old

print "\nThe example of calculating of Coulomb interactions (the method of Ewald summation)"

# reading the particle coordinates, charges and box size from old espresso data file
# file 'Coulomb_Ewald.dat' contains the data we need (Lx, Ly, Lz, x, y, z, type, q)
print "Reading system data:"
Lx, Ly, Lz, x, y, z, type, q, vx,vy,vz,fx,fy,fz,bondpairs = espresso_old.read('Coulomb_Ewald.dat')

# creating the system box
box = (Lx, Ly, Lz)
print "System box size:", box
# number of particles
num_particles = len(x)
print "Number of particles =  ", num_particles

'''
#  Ewald method suppose to calculate electrostatic interaction dividing it into R space and
#  K space part
#  
#  alpha - Ewald parameter
#  rspacecutoff - the cutoff in real space
#  kspacecutoff - the cutoff in reciprocal space   
'''
alpha          = 1.2
print "Ewald parameter:       ", alpha
rspacecutoff   = 2.5
print "The cutoff in R space: ", rspacecutoff
kspacecutoff   = 30
print "The cutoff in K space: ", kspacecutoff

# seting the skin for Verlet list
skin           = 0.3

# Coulomb prefactor parameter is the product of Bjerrum length and temperature
bjerrumlength  = 1.0
temperature    = 1.0
coulomb_prefactor = bjerrumlength * temperature
print "Coulomb prefactor:     ", coulomb_prefactor

nodeGrid       = espresso.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, rspacecutoff, skin)
system         = espresso.System()
system.rng     = espresso.esutil.RNG()
system.bc      = espresso.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# Adding the particles
props = ['id', 'pos', 'type', 'q']
new_particles = []
for i in range(0, num_particles):
  part = [ i, Real3D(x[i], y[i], z[i]), type[i], q[i] ]
  new_particles.append(part)
system.storage.addParticles(new_particles, *props)
system.storage.decompose()

## potential and interaction ##

# the R space part of electrostatic interaction according to the Ewald method
# setting the Verlet list
vl = espresso.VerletList(system, rspacecutoff)
'''
  Creating the Coulomb potential which is responsible for the R space part according to the
  Ewald method.
  It is based on the Coulomb prefactor (coulomb_prefactor), Ewald parameter (alpha),
  and the cutoff in R space (rspacecutoff)
'''
coulombR_pot = espresso.interaction.CoulombRSpace(coulomb_prefactor, alpha, rspacecutoff)
# creating the interaction based on the Verlet list
coulombR_int = espresso.interaction.VerletListCoulombRSpace(vl)
# setting the potential for the interaction between particles of type 0 and 0
coulombR_pot = coulombR_int.setPotential(type1=0, type2=0, potential = coulombR_pot)
# adding the interaction to the system
system.addInteraction(coulombR_int)

# the K space part of electrostatic interaction according to the Ewald method
'''
  Creating the Coulomb potential which is responsible for the K space part according to the
  Ewald method.
  It is based on the system information (system), Coulomb prefactor (coulomb_prefactor),
  Ewald parameter (alpha), and the cutoff in K space (kspacecutoff)
'''
ewaldK_pot = espresso.interaction.CoulombKSpaceEwald(system, coulomb_prefactor, alpha, kspacecutoff)
# creating the interaction based on the Cell list for all particle interaction and potential in K space
ewaldK_int = espresso.interaction.CellListCoulombKSpaceEwald(system.storage, ewaldK_pot)
# adding the interaction to the system
system.addInteraction(ewaldK_int)

# creating the integrator
integrator    = espresso.integrator.VelocityVerlet(system)
# seting the time step
integrator.dt = 0.005

# analysis
temperature = espresso.analysis.Temperature(system)

pr_format = '%6d %8.4f %10.8f %10.8f %10.8f  %10.8f'
T = temperature.compute()
Ek = 0.5 * T * (3 * num_particles)
Ersp = coulombR_int.computeEnergy()
Eksp = ewaldK_int.computeEnergy()
Etotal = Ek + Ersp + Eksp
print('\n  step      T        Etotal      Ek        Ersp         Eksp\n')
print(pr_format % (0, T, Etotal, Ek, Ersp, Eksp))

nsteps = 5
intervals = 10

for i in range (1, intervals+1):
  integrator.run(nsteps)
  step = nsteps * i
  T = temperature.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ersp = coulombR_int.computeEnergy()
  Eksp = ewaldK_int.computeEnergy()
  Etotal = Ek + Ersp + Eksp
  print(pr_format % (step, T, Etotal, Ek, Ersp, Eksp))
