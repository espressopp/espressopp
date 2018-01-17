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

#########################################################################################
#                                                                                       #
#  ESPResSo++ Python script for an H-AdResS tetrahedral liquid including FEC pressure   #
#                                                                                       #
#  FEC stands for Free Energy Correction                                                #
#########################################################################################

import sys
import time
import espressopp
import mpi4py.MPI as MPI

import Tetracryst # preparation of tetrahedral crystal and constuctions of bonds in tetrahedral liquid

from espressopp import Real3D, Int3D
from espressopp.tools import decomp
from espressopp.tools import timers

# integration steps, cutoff, skin, AdResS specifications
steps = 1000
timestep = 0.0005
intervals = 100

rc = 4.5 # cutoff coarse-grained potential
rca = 1.122462048309373 # cutoff atomistic potential (cutoff (2^(1/6)), WCA)
skin = 0.4

# parameters for the thermostat
#gamma = 2.0
#temp = 1.0

# parameters for size of AdResS dimensions
ex_size = 5.0
hy_size = 5.0

# read equilibrated configuration file
pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz("equilibrated_conf.xyz")

# Table for coarse-grained potential
tabCG = "table_potential.dat"

# FEC compensation table
tabFEC = "table_FEC_Helmholtz.dat"

# number of CG particles
num_particlesCG = len(x)/4

# number of AT particles
num_particles = len(x)

# set up the system
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

# (H-)AdResS domain decomposition
system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)


# prepare AT particles
allParticlesAT = []
allParticles = []
tuples = []
for pidAT in range(num_particles):
    allParticlesAT.append([pidAT, # add here these particles just temporarily
                         Real3D(x[pidAT], y[pidAT], z[pidAT]), # position
                         Real3D(vx[pidAT], vy[pidAT], vz[pidAT]), # velocity
                         Real3D(0, 0, 0), # force
                         1, 1.0, 1]) # type, mass, is AT particle

# create CG particles
for pidCG in range(num_particlesCG):
    # we put CG molecule in first atom, later CG molecules will be positioned in the center
    cmp = espressopp.tools.AdressSetCG(4, pidCG, allParticlesAT)
    # Preparation of tuples (tuples define, which atoms belong to which CG molecules)
    tmptuple = [pidCG+num_particles]
    for pidAT2 in range(4):
        pid = pidCG*4+pidAT2
        tmptuple.append(pid)

    # append CG particles
    allParticles.append([pidCG+num_particles, # CG particle has to be added first!
                         Real3D(cmp[0], cmp[1], cmp[2]), # pos
                         Real3D(0, 0, 0), # vel
                         Real3D(0, 0, 0), # force
                         0, 4.0, 0]) # type, mass, is not AT particle
    # append AT particles
    for pidAT in range(4):
        pid = pidCG*4+pidAT
        allParticles.append([pid, # now the AT particles can be added
                            (allParticlesAT[pid])[1], # pos
                            (allParticlesAT[pid])[2], # vel
                            (allParticlesAT[pid])[3], # force
                            (allParticlesAT[pid])[4], # type
                            (allParticlesAT[pid])[5], # mass
                            (allParticlesAT[pid])[6]]) # is AT particle
    # append tuple to tuplelist
    tuples.append(tmptuple)


# add particles to system
system.storage.addParticles(allParticles, "id", "pos", "v", "f", "type", "mass", "adrat")

# create FixedTupleList object
ftpl = espressopp.FixedTupleListAdress(system.storage)

# and add the tuples
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)

# add bonds between AT particles
fpl = espressopp.FixedPairListAdress(system.storage, ftpl)
bonds = Tetracryst.makebonds(len(x))
fpl.addBonds(bonds)

# decompose after adding tuples and bonds
print "Added tuples and bonds, decomposing now ..."
system.storage.decompose()
print "done decomposing"

# AdResS Verlet list
vl = espressopp.VerletListAdress(system, cutoff=rc, adrcut=rc,
                                dEx=ex_size, dHy=hy_size,
                                adrCenter=[Lx/2, Ly/2, Lz/2])

# non-bonded potentials
# LJ capped WCA between AT and tabulated potential between CG particles
interNB = espressopp.interaction.VerletListHadressLennardJones(vl, ftpl) # Here we need specific (H-)AdResS interaction type
potWCA  = espressopp.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift='auto', cutoff=rca)
potCG = espressopp.interaction.Tabulated(itype=3, filename=tabCG, cutoff=rc) # CG
interNB.setPotentialAT(type1=1, type2=1, potential=potWCA) # AT
interNB.setPotentialCG(type1=0, type2=0, potential=potCG) # CG
system.addInteraction(interNB)

# bonded potentials
# quartic potential between AT particles
potQuartic = espressopp.interaction.Quartic(K=75.0, r0=1.0)
interQuartic = espressopp.interaction.FixedPairListQuartic(system, fpl, potQuartic)
system.addInteraction(interQuartic)

# VelocityVerlet integrator
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = timestep

# add AdResS extension
adress = espressopp.integrator.Adress(system, vl, ftpl)
integrator.addExtension(adress)

# add Langevin thermostat extension
#langevin = espressopp.integrator.LangevinThermostat(system)
#langevin.gamma = gamma
#langevin.temperature = temp
#langevin.adress = True # enable AdResS!
#integrator.addExtension(langevin)

# add FEC
fec = espressopp.integrator.FreeEnergyCompensation(system, center=[Lx/2, Ly/2, Lz/2])
fec.addForce(itype=3, filename=tabFEC, type=0)
integrator.addExtension(fec)

# distribute atoms and CG molecules according to AdResS domain decomposition, place CG molecules in the center of mass
espressopp.tools.AdressDecomp(system, integrator)

# system information
print ''
print 'AdResS Center =', [Lx/2, Ly/2, Lz/2]
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

fmt = '%5d %8.4f %12.3f %12.3f %12.3f %12.3f %12.3f\n'
T = temperature.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interNB.computeEnergy()
Eb = interQuartic.computeEnergy()
Ecorr = fec.computeCompEnergy()
sys.stdout.write(' step    Temp       etotal      enonbonded    ebonded     ekinetic     ecorrection\n')
sys.stdout.write(fmt % (0, T, Ek + Ep + Eb + Ecorr, Ep, Eb, Ek, Ecorr))

# pressure profile preparation
pressure_array_total = []
Adds = 0.0
pressureprofilegrid = 100
pressureprofile = espressopp.analysis.XPressure(system)

# timer, steps
nsteps = steps / intervals
start_time = time.clock()

# integration and on the fly analysis
for s in range(1, intervals + 1):
  integrator.run(nsteps)
  step = nsteps * s
  T = temperature.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interNB.computeEnergy()
  Eb = interQuartic.computeEnergy()
  Ecorr = fec.computeCompEnergy()
  sys.stdout.write(fmt % (step, T, Ek + Ep + Eb + Ecorr, Ep, Eb, Ek, Ecorr))

  # calculate pressure profile
  pressure_array = pressureprofile.compute(pressureprofilegrid)
  for i in range(len(pressure_array)):
    if(i>=len(pressure_array_total)):
      pressure_array_total.append(pressure_array[i])
    else:
      pressure_array_total[i] += pressure_array[i]
  Adds += 1.0

# correct the pressure profile according to number of samples
for i in range(len(pressure_array_total)):
  pressure_array_total[i] /= Adds

# printing pressure profile
nameFile = 'pressure_profile_Helmholtz.dat'
print ''
print "Printing the pressure profile to %s\n" %nameFile
tempFile = open (nameFile, 'w')
fmt = ' %12.8f %12.8f\n'
dr = Lx / float(pressureprofilegrid)
for i in range( len(pressure_array_total) ):
  tempFile.write(fmt % ( (i+0.5)*dr, pressure_array_total[i] ))
tempFile.close()

# simulation information
end_time = time.clock()
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))
