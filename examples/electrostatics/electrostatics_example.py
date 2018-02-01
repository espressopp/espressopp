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

#############################################################################
#                                                                           #
#  ESPResSo++ Python script for the calculation of long range interactions  #
#  (Coulomb interaction) using the Ewald summation and the P3M methods.     #
#                                                                           #
#############################################################################
'''
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

import mpi4py.MPI as MPI
import espressopp
from espressopp import Real3D

# initial parameters
N = 16                 # number of particles on lattice site
num_particles = N**3   # total number of particles 
rho = 0.03              # number density of particles, number of particles devided by volume

# creating a cubic NaCl crystal
#print 'Creating a simple cubic structure...'
x, y, z, Lx, Ly, Lz = espressopp.tools.lattice.createCubic(num_particles, rho, False)

# creating the system box
box = (Lx, Ly, Lz)
print 'System box size: ',      box
print 'Number of particles = ', num_particles

#  Ewald summation parameters
#alphaEwald     = 1.112583061 #  alpha - Ewald parameter
alphaEwald     = 0.660557
rspacecutoff   = 4.9 #  rspacecutoff - the cutoff in real space
kspacecutoff   = 30 #  kspacecutoff - the cutoff in reciprocal space

print 'Ewald parameters:'
print 'alfa=%f, rcutoff=%f, kcutoff=%d' % (alphaEwald, rspacecutoff, kspacecutoff)

#  P3M parameters
M              = espressopp.Int3D(64, 64, 64)
P              = 7
#alphaP3M       = 1.112583061 #  alpha - Ewald parameter
alphaP3M       = 0.660557

print 'P3M parameters:'
print 'Mesh=', M,', charge assignment order=%d,  alphaP3M=%lf' % ( P, alphaP3M)

# a skin for Verlet list
skin              = 0.2

# Coulomb prefactor parameters
bjerrumlength     = 1.0
temperature       = 1.0
coulomb_prefactor = bjerrumlength * temperature

nodeGrid          = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size,box,rspacecutoff,skin)
cellGrid          = espressopp.tools.decomp.cellGrid(box, nodeGrid, rspacecutoff, skin)

print ''
print 'density = %.4f' % (rho)
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

'''
  Below two systems for Ewald summation and PPPM methods will be created.
'''

#######################################################################################
#   system for Ewald 
#######################################################################################

systemEwald         = espressopp.System()
systemEwald.rng     = espressopp.esutil.RNG()
systemEwald.bc      = espressopp.bc.OrthorhombicBC(systemEwald.rng, box)
systemEwald.skin    = skin
systemEwald.storage = espressopp.storage.DomainDecomposition(systemEwald, nodeGrid, cellGrid)

#######################################################################################
#   system for PPPM
#######################################################################################

systemPPPM         = espressopp.System()
systemPPPM.rng     = espressopp.esutil.RNG()
systemPPPM.bc      = espressopp.bc.OrthorhombicBC(systemPPPM.rng, box)
systemPPPM.skin    = skin
systemPPPM.storage = espressopp.storage.DomainDecomposition(systemPPPM, nodeGrid, cellGrid)

#######################################################################################

# adding particles
props = ['id', 'pos', 'type', 'q']
new_particles = []
countX = countY = countZ = 0
for i in range(0, num_particles):
  
  # charge should be accordingly to NaCl crystall
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
vlEwald = espressopp.VerletList(systemEwald, rspacecutoff)
vlPPPM  = espressopp.VerletList(systemPPPM,  rspacecutoff)

# real space interaction for Ewald system
# R space part of electrostatic interaction
coulombR_potEwald = espressopp.interaction.CoulombRSpace(coulomb_prefactor, alphaEwald, rspacecutoff)
# creating an interaction based on the Verlet list
coulombR_intEwald = espressopp.interaction.VerletListCoulombRSpace(vlEwald)
# setting the potential for the interaction between particles of type 0 and 0
coulombR_intEwald.setPotential(type1=0, type2=0, potential = coulombR_potEwald)
# adding the interaction to the system
systemEwald.addInteraction(coulombR_intEwald)

# real space interaction for PPPM system
# R space part of electrostatic interaction
coulombR_potP3M = espressopp.interaction.CoulombRSpace(coulomb_prefactor, alphaP3M, rspacecutoff)
# creating an interaction based on the Verlet list
coulombR_intPPPM = espressopp.interaction.VerletListCoulombRSpace(vlPPPM)
# setting the potential for the interaction between particles of type 0 and 0
coulombR_intPPPM.setPotential(type1=0, type2=0, potential = coulombR_potP3M)
# adding the interaction to the system
systemPPPM.addInteraction(coulombR_intPPPM)

# K space part of electrostatic interaction
ewaldK_pot = espressopp.interaction.CoulombKSpaceEwald(systemEwald, coulomb_prefactor, alphaEwald, kspacecutoff)
# creating an interaction based on the Cell list for all particle interaction and potential in K space
ewaldK_int = espressopp.interaction.CellListCoulombKSpaceEwald(systemEwald.storage, ewaldK_pot)
# adding the interaction to the system
systemEwald.addInteraction(ewaldK_int)

# PPPM system
p3m_pot = espressopp.interaction.CoulombKSpaceP3M( systemPPPM, coulomb_prefactor, alphaP3M, M, P, rspacecutoff)
# creating the interaction based on the Cell list for all particle interaction and potential in K space
p3m_int = espressopp.interaction.CellListCoulombKSpaceP3M(systemPPPM.storage, p3m_pot)
# adding the interaction to the system
systemPPPM.addInteraction(p3m_int)

### Integrators for Ewald and PPPM
# creating the integrator which based on the Verlet algorithm
integratorEwald    = espressopp.integrator.VelocityVerlet(systemEwald)
# seting the time step (it is not important here)
integratorEwald.dt = 0.0001
# nothing will be changed in system, just forces will be calculated ones
integratorEwald.run(0)

# creating the integrator which based on the Verlet algorithm
integratorPPPM    = espressopp.integrator.VelocityVerlet(systemPPPM)
# seting the time step (it is not important here)
integratorPPPM.dt = 0.0001
# nothing will be changed in system, just forces will be calculated ones
integratorPPPM.run(0)

# printing the particle id and force difference (x,y,z) for first 6 particles
print ('\n    Difference between forces calculated by Ewald summation and PPPM (first 6 particles)')
print ('%3s %20s %20s %20s\n' % ('id', 'dfx', 'dfy', 'dfz'))

#sock = espressopp.tools.vmd.connect(systemPPPM)
#espressopp.tools.vmd.imd_positions(systemPPPM, sock)

print_N = min(num_particles, 20)

for j in range(0, print_N):
  print ( '%3d     %3.17f     %3.17f     %3.17f' % (j, \
    abs(systemEwald.storage.getParticle(j).f.x - systemPPPM.storage.getParticle(j).f.x), \
    abs(systemEwald.storage.getParticle(j).f.y - systemPPPM.storage.getParticle(j).f.y), \
    abs(systemEwald.storage.getParticle(j).f.z - systemPPPM.storage.getParticle(j).f.z)) )
  
  print 'force:', systemPPPM.storage.getParticle(j).f, '   ', systemEwald.storage.getParticle(j).f
  

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
