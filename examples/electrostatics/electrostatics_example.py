'''
#  This script is an example of calculation of long range interactions
#  (Coulomb interaction) using the Ewald summation and the P3M methods.
#  
#  Initially, the simple cubic structure is generated in order to represent the
#  NaCl crystal. Then the energy and forces are calculated and compared using both
#  the Ewald summation and the P3M. At the end the Madelung constant of NaCl crystal
#  is calculated.
#  
#  At the moment there is only metallic surrounding media is possible.

#  Parameters:

#  Ewald summation:
    alpha          = 1.112583061   (Ewald parameter)
    rspacecutoff   = 4.9           (the cutoff in real space)
    kspacecutoff   = 30            (the cutoff in reciprocal space)

#  P3M:
    M              = (16, 16, 16)  (mesh)
    P              = 7             (charge assignment order)
'''

# The script itself
import MPI
import espresso
from espresso import Real3D

# initial parameters
N = 2                 # number of particles on lattice site
num_particles = N**3   # total number of particles 
rho = 0.0079999999              # number density of particles, number of particles devided by volume

# creating a cubic NaCl crystal
#print 'Creating a simple cubic structure...'
#x, y, z, Lx, Ly, Lz = espresso.tools.init_cfg.lattice.createCubic(num_particles, rho)

x = []
y = []
z = []
for i in range(0,8):
  x.append(0.0)
  y.append(0.0)
  z.append(0.0)

x[0]=4.159994;   y[0]=0.919649;   z[0]=7.564105
x[1]=5.297002;   y[1]=9.304365;   z[1]=3.835021
x[2]=6.539190;   y[2]=0.668422;   z[2]=7.226604
x[3]=6.711494;   y[3]=3.834157;   z[3]=6.316347
x[4]=8.847071;   y[4]=5.194164;   z[4]=6.515186
x[5]=2.377744;   y[5]=2.624530;   z[5]=7.621980
x[6]=7.533558;   y[6]=9.092081;   z[6]=0.726859
x[7]=2.727100;   y[7]=8.976563;   z[7]=2.749068

Lx = Ly = Lz = 10.0

# creating the system box
box = (Lx, Ly, Lz)
print 'System box size: ',      box
print 'Number of particles = ', num_particles

#  Ewald summation parameters
alpha          = 1.112583061 #  alpha - Ewald parameter
rspacecutoff   = 4.9 #  rspacecutoff - the cutoff in real space
kspacecutoff   = 30 #  kspacecutoff - the cutoff in reciprocal space

print 'Ewald parameters:'
print 'alfa=%f, rcutoff=%f, kcutoff=%d' % (alpha, rspacecutoff, kspacecutoff)

#  P3M parameters
M              = espresso.Int3D(16, 16, 16)
P              = 7

print 'P3M parameters:'
print 'Mesh=', M, ', charge assignment order=%d' % (P)

# a skin for Verlet list
skin           = 0.2

# Coulomb prefactor parameters
bjerrumlength     = 1.0
temperature       = 1.0
coulomb_prefactor = bjerrumlength * temperature

nodeGrid          = espresso.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
cellGrid          = espresso.tools.decomp.cellGrid(box, nodeGrid, rspacecutoff, skin)

'''
  Below two systems for Ewald summation and PPPM methods will be created.
'''

#######################################################################################
#   system for Ewald 
#######################################################################################

systemEwald         = espresso.System()
systemEwald.rng     = espresso.esutil.RNG()
systemEwald.bc      = espresso.bc.OrthorhombicBC(systemEwald.rng, box)
systemEwald.skin    = skin
systemEwald.storage = espresso.storage.DomainDecomposition(systemEwald, nodeGrid, cellGrid)

#######################################################################################
#   system for PPPM
#######################################################################################

systemPPPM         = espresso.System()
systemPPPM.rng     = espresso.esutil.RNG()
systemPPPM.bc      = espresso.bc.OrthorhombicBC(systemPPPM.rng, box)
systemPPPM.skin    = skin
systemPPPM.storage = espresso.storage.DomainDecomposition(systemPPPM, nodeGrid, cellGrid)

#######################################################################################

# adding particles
props = ['id', 'pos', 'type', 'q']
new_particles = []
countX = countY = countZ = 0
for i in range(0, num_particles):
  
  charge = pow(-1, countX + countY + countZ)
  part = [ i, Real3D(x[i], y[i], z[i]), 0, charge ]
  new_particles.append(part)
  
  countX += 1
  if countX >= N:
    countX = 0
    countY += 1
    if countY >= N:
      countY = 0
      countZ += 1


# adding particles to Ewald system
systemEwald.storage.addParticles(new_particles, *props)
systemEwald.storage.decompose()

# adding particles to PPPM system
systemPPPM.storage.addParticles(new_particles, *props)
systemPPPM.storage.decompose()

## potentials and interactions ##

# setting a Verlet list
vlEwald = espresso.VerletList(systemEwald, rspacecutoff)
vlPPPM = espresso.VerletList(systemPPPM, rspacecutoff)

# R space part of electrostatic interaction
# potential is the same for Ewald summation and PPPM
coulombR_pot = espresso.interaction.CoulombRSpace(coulomb_prefactor, alpha, rspacecutoff)

# real space interaction for Ewald system
# creating an interaction based on the Verlet list
coulombR_intEwald = espresso.interaction.VerletListCoulombRSpace(vlEwald)
# setting the potential for the interaction between particles of type 0 and 0
coulombR_intEwald.setPotential(type1=0, type2=0, potential = coulombR_pot)
# adding the interaction to the system
systemEwald.addInteraction(coulombR_intEwald)

# real space interaction for PPPM system
# creating an interaction based on the Verlet list
coulombR_intPPPM = espresso.interaction.VerletListCoulombRSpace(vlPPPM)
# setting the potential for the interaction between particles of type 0 and 0
coulombR_intPPPM.setPotential(type1=0, type2=0, potential = coulombR_pot)
# adding the interaction to the system
systemPPPM.addInteraction(coulombR_intPPPM)

# K space part of electrostatic interaction
ewaldK_pot = espresso.interaction.CoulombKSpaceEwald(systemEwald, coulomb_prefactor, alpha, kspacecutoff)
# creating an interaction based on the Cell list for all particle interaction and potential in K space
ewaldK_int = espresso.interaction.CellListCoulombKSpaceEwald(systemEwald.storage, ewaldK_pot)
# adding the interaction to the system
systemEwald.addInteraction(ewaldK_int)

# PPPM system
p3m_pot = espresso.interaction.CoulombKSpaceP3M( systemPPPM, coulomb_prefactor, alpha, M, P, rspacecutoff)
# creating the interaction based on the Cell list for all particle interaction and potential in K space
p3m_int = espresso.interaction.CellListCoulombKSpaceP3M(systemPPPM.storage, p3m_pot)
# adding the interaction to the system
systemPPPM.addInteraction(p3m_int)

hhh = ( p3m_int.computeEnergy() + coulombR_intPPPM.computeEnergy() )
### Integrators for Ewald and PPPM
print '   PPP_energy: ', hhh, p3m_int.computeEnergy(), coulombR_intPPPM.computeEnergy()
print '   Ewald_energy: ', ewaldK_int.computeEnergy(), coulombR_intEwald.computeEnergy()

# creating the integrator which based on the Verlet algorithm
integratorEwald    = espresso.integrator.VelocityVerlet(systemEwald)
# seting the time step (it is not important here)
integratorEwald.dt = 0.0001
# nothing will be changed in system, just forces will be calculated ones
integratorEwald.run(0)

# creating the integrator which based on the Verlet algorithm
integratorPPPM    = espresso.integrator.VelocityVerlet(systemPPPM)
# seting the time step (it is not important here)
integratorPPPM.dt = 0.0001
# nothing will be changed in system, just forces will be calculated ones
integratorPPPM.run(0)

# printing the particle id and force difference (x,y,z) for first 6 particles
print ('\n    Difference between forces calculated by Ewald summation and PPPM (first 6 particles)')
print ('%3s %20s %20s %20s\n' % ('id', 'dfx', 'dfy', 'dfz'))

#sock = espresso.tools.vmd.connect(systemPPPM)
#espresso.tools.vmd.imd_positions(systemPPPM, sock)

for j in range(0, num_particles):
  print 'position:', j, systemPPPM.storage.getParticle(j).pos, '   ', systemEwald.storage.getParticle(j).pos
  
  
for j in range(0, num_particles):
  print ( '%3d     %3.17f     %3.17f     %3.17f' % (j, \
    abs(systemEwald.storage.getParticle(j).f.x - systemPPPM.storage.getParticle(j).f.x), \
    abs(systemEwald.storage.getParticle(j).f.y - systemPPPM.storage.getParticle(j).f.y), \
    abs(systemEwald.storage.getParticle(j).f.z - systemPPPM.storage.getParticle(j).f.z)) )
  
  print 'force:', systemPPPM.storage.getParticle(j).f, '   ', systemEwald.storage.getParticle(j).f

#print 'no energy calc'
#exit(0)

# calculating the R space part of electrostatic energy
energyEwaldR = coulombR_intEwald.computeEnergy()
# calculating the K space part of electrostatic energy
energyEwaldK = ewaldK_int.computeEnergy()
# total energy (Ewald summation)
enTotEwald = energyEwaldR + energyEwaldK


# calculating the R space part of electrostatic energy
energyPPPMR = coulombR_intPPPM.computeEnergy()
# calculating the K space part of electrostatic energy
energyPPPMK = p3m_int.computeEnergy()
# total energy (PPPM)
enTotPPPM = energyPPPMR + energyPPPMK



# printing the total energy and the difference
print 'Energy (Ewald summation): %5.16f  Energy (PPPM): %5.16f\n' % (enTotEwald, enTotPPPM)
print 'The difference in energy (Ewald - PPPM): %5.16f\n' % (enTotEwald-enTotPPPM)

a = 2 * pow( Lx*Ly*Lz / num_particles , 1./3. )
madelung_NaCl = -1.747564594633182190636212035544397403481
print ("Madelung constant is: %14.10f\n" % (enTotEwald/num_particles * a))
print (" error: %e\n\n" % ( abs( abs( enTotPPPM/num_particles * a) - abs(madelung_NaCl))))
