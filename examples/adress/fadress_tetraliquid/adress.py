#!/usr/bin/env python2 
#  Copyright (C) 2016, 2017(H)
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
#  ESPResSo++ Python script for an F-AdResS tetrahedral liquid simulation #
#                                                                         #
###########################################################################

import sys
import time
import espressopp
import mpi4py.MPI as MPI
import logging
from espressopp import Real3D, Int3D
from espressopp.tools import espresso_old
from espressopp.tools import decomp
from espressopp.tools import timers

# integration steps, cutoff, skin and thermostat flag (nvt = False is nve)
steps = 5000
timestep = 0.0001
intervals = 500

rc = 2.31 # CG cutoff, Morse
rca = 1.122462048309373 # AT cutoff (2^(1/6)), WCA
skin = 0.4

gamma = 0.5
temp = 1.0

ex_size = 12.5
hy_size = 5.0

# writes the tabulated file
def writeTabFile(pot, name, N, low=0.0, high=2.5, body=2):
    outfile = open(name, "w")
    delta = (high - low) / (N - 1)

    for i in range(N):
        r = low + i * delta
        energy = pot.computeEnergy(r)
        if body == 2:# this is for 2-body potentials
            force = pot.computeForce(Real3D(r, 0.0, 0.0))[0]
            #force /= r
        else: # this is for 3- and 4-body potentials
            force = pot.computeForce(r)
        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))

    outfile.close()

# tabulated morse potential used for CG interactions
tabMorse = "pot-morse.txt"
potMorse = espressopp.interaction.Morse(epsilon=0.105, alpha=2.4, rMin=rc, cutoff=rc, shift="auto")
writeTabFile(potMorse, tabMorse, N=512, low=0.005, high=4.5)

# read ESPResSo configuration file
Lx, Ly, Lz, x, y, z, type, q, vx, vy, vz, fx, fy, fz, bonds = espresso_old.read("adress.espressopp")
num_particlesCG = 5001 # number of VP/CG particles
num_particles = len(x) # 20004

sys.stdout.write('Setting up simulation ...\n')
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

system = espressopp.System()
system.rng = espressopp.esutil.RNG()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size,size,rc,skin)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)

# AdResS domain decomposition
system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)


# prepare AT particles
allParticlesAT = []
allParticles = []
tuples = []
for pidAT in range(num_particles):
    allParticlesAT.append([pidAT, # add here these particles just temporarly!
                         Real3D(x[pidAT], y[pidAT], z[pidAT]),
                         Real3D(vx[pidAT], vy[pidAT], vz[pidAT]),
                         Real3D(fx[pidAT], fy[pidAT], fz[pidAT]),
                         1, 1.0, 1]) # type, mass, is AT particle

# create CG particles from center of mass
for pidCG in range(num_particlesCG):
    cmp = [0,0,0]
    cmv = [0,0,0]
    tmptuple = [pidCG+num_particles]
    # com calculation
    for pidAT in range(4):
        pid = pidCG*4+pidAT
        tmptuple.append(pid)
        pos = (allParticlesAT[pid])[1]
        vel = (allParticlesAT[pid])[2]
        for i in range(3):
            cmp[i] += pos[i] # masses are 1.0 so we skip multiplication
            cmv[i] += vel[i]
    for i in range(3):
        cmp[i] /= 4.0 # 4.0 is the mass of molecule
        cmv[i] /= 4.0

    allParticles.append([pidCG+num_particles, # CG particle has to bo added first!
                         Real3D(cmp[0], cmp[1], cmp[2]), # pos
                         Real3D(cmv[0], cmv[1], cmv[2]), # vel
                         Real3D(0, 0, 0), # f
                         0, 4.0, 0]) # type, mass, is not AT particle

    for pidAT in range(4):
        pid = pidCG*4+pidAT
        allParticles.append([pid, # now the AT particles can be added
                            (allParticlesAT[pid])[1], # pos
                            (allParticlesAT[pid])[2], # vel
                            (allParticlesAT[pid])[3], # f
                            (allParticlesAT[pid])[4], # type
                            (allParticlesAT[pid])[5], # mass
                            (allParticlesAT[pid])[6]]) # is AT particle

    tuples.append(tmptuple)

# add particles
system.storage.addParticles(allParticles, "id", "pos", "v", "f", "type", "mass", "adrat")

# add tuples
ftpl = espressopp.FixedTupleListAdress(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)

# add bonds between AT particles
fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
fpl.addBonds(bonds)

# decompose after adding tuples and bonds
print "Added tuples and bonds, decomposing now ..."
system.storage.decompose()

print "done decomposing"

# AdResS Verlet list
vl = espressopp.VerletListAdress(system, cutoff=rc+skin, adrcut=rc+skin,
                                dEx=ex_size, dHy=hy_size,
                                adrCenter=[18.42225, 18.42225, 18.42225])

# non-bonded potentials
# LJ Capped WCA between AT and tabulated Morse between CG particles
interNB = espressopp.interaction.VerletListAdressLennardJonesCapped(vl, ftpl)
potWCA  = espressopp.interaction.LennardJonesCapped(epsilon=1.0, sigma=1.0, shift=True, caprad=0.27, cutoff=rca)
potMorse = espressopp.interaction.Tabulated(itype=2, filename=tabMorse, cutoff=rc) # CG
interNB.setPotentialAT(type1=1, type2=1, potential=potWCA) # AT
interNB.setPotentialCG(type1=0, type2=0, potential=potMorse) # CG
system.addInteraction(interNB)

# bonded potentials
# FENE and LJ potential between AT particles
potFENE = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
potLJ = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift=True, cutoff=rca)
interFENE = espressopp.interaction.FixedPairListFENE(system, fpl, potFENE)
interLJ = espressopp.interaction.FixedPairListLennardJones(system, fpl, potLJ)
system.addInteraction(interFENE)
system.addInteraction(interLJ)

# VV integrator
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = timestep

# add AdResS extension
adress = espressopp.integrator.Adress(system, vl,ftpl)
integrator.addExtension(adress)

# add Langevin thermostat extension
langevin = espressopp.integrator.LangevinThermostat(system)
langevin.gamma = gamma
langevin.temperature = temp
langevin.adress = True # enable AdResS!
integrator.addExtension(langevin)


print ''
print 'number of AT particles =', num_particles
print 'number of CG particles =', num_particlesCG
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'steps =', steps
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
temperature = espressopp.analysis.Temperature(system)
pressure = espressopp.analysis.Pressure(system)
pressureTensor = espressopp.analysis.PressureTensor(system)

fmt = '%5d %8.4f %10.5f %8.5f %12.3f %12.3f %12.3f %12.3f\n'

T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interNB.computeEnergy()
Eb = interFENE.computeEnergy()
sys.stdout.write(' step     T          P        Pxy       etotal     epotential      ebonded     ekinetic\n')
sys.stdout.write(fmt % (0, T, P, Pij[3], Ek + Ep + Eb, Ep, Eb, Ek))

start_time = time.clock()
nsteps = steps / intervals
for s in range(1, intervals + 1):
  integrator.run(nsteps)
  step = nsteps * s
  T = temperature.compute()
  P = pressure.compute()
  Pij = pressureTensor.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interNB.computeEnergy()
  Eb = interFENE.computeEnergy()
  sys.stdout.write(fmt % (step, T, P, Pij[3], Ek + Ep + Eb, Ep, Eb, Ek))
  system.storage.decompose()

end_time = time.clock()

# simulation information
end_time = time.clock()
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))

