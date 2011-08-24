#!/usr/bin/env python
# -*- coding: utf-8 -*-



###########################################################################
#                                                                         #
#  This Python script may be used to simulate a monatomic LJ fluid in the #
#  NVE or NVT ensemble. The starting configuration may be taken from      #
#  either a LAMMPS data file or by generating coordinates on a lattice.   #
#                                                                         #
###########################################################################

import sys
import time
import espresso
import MPI
import logging
from espresso import Real3D, Int3D
from espresso.tools.convert import espresso_old
from espresso.tools import decomp
from espresso.tools import timers

# integration steps, cutoff, skin and thermostat flag (nvt = False is nve)
steps = 100
rc = 2.31 # CG cutoff, Morse
rca = 1.122462048309373 # AT cutoff (2^(1/6)), WCA
skin = 0.8
nvt = True
timestep = 0.01
intervals = 10




# writes the tabulated file
def writeTabFile(pot, name, N, low=0.0, high=2.5, body=2):
    outfile = open(name, "w")
    delta = (high - low) / (N - 1)
     
    for i in range(N):
        r = low + i * delta
        energy = pot.computeEnergy(r)
        if body == 2:# this is for 2-body potentials
            force = pot.computeForce(Real3D(r, 0.0, 0.0))[0]
        else: # this is for 3- and 4-body potentials
            force = pot.computeForce(r)
        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))
     
    outfile.close()



tabWCA = "pot-wca.txt"
potLJ  = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift=True, cutoff=rca)
writeTabFile(potLJ, tabWCA, N=512, low=0.005, high=rca)

#tabFENE = "pot-fene.txt"
#potFENE  = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
#writeTabFile(potFENE, tabFENE, N=512, low=0.005, high=2.0)

tabMorse = "pot-morse.txt"
potMorse = espresso.interaction.Morse(epsilon=0.105, alpha=2.4, rMin=rc, cutoff=rc, shift="auto")
writeTabFile(potMorse, tabMorse, N=512, low=0.005, high=rc)


# read ESPResSo configuration file 
Lx, Ly, Lz, x, y, z, type, q, vx, vy, vz, fx, fy, fz, bonds = espresso_old.read("ex_j1.txt")
num_particlesCG = 5000 # number of VP/CG particles
num_particles = len(x) - num_particlesCG  # 20004 = 25004 - 5000 


######################################################################
### IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE ###
######################################################################



sys.stdout.write('Setting up simulation ...\n')
print Lx, Ly, Lz
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)
system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

# compute the number of cells on each node
def calcNumberCells(size, nodes, cutoff):
    ncells = 1
    while size / (ncells * nodes) >= cutoff:
       ncells = ncells + 1
    return ncells - 1

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)

system.storage = espresso.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)


# add CG and AT particles
allParticles = []
tuples = []
idx = 0 # just an index
for pidCG in range(num_particles-1, num_particles+num_particlesCG):
    #print pidCG,
    tmptuple = [pidCG] # first pid in tuple is pid of CG/VP particle
    allParticles.append([pidCG,
                         Real3D(x[pidCG], y[pidCG], z[pidCG]),
                         Real3D(vx[pidCG], vy[pidCG], vz[pidCG]),
                         Real3D(fx[pidCG], fy[pidCG], fz[pidCG]), 0, 4, 0])
    
    for pidAT in range(4): # each CG has 4 AT particles
        #print pidAT+idx,
        tmptuple.append(pidAT+idx) # other pids in tuple are pids of AT particles
        allParticles.append([pidAT+idx, Real3D(x[pidAT+idx], y[pidAT+idx], z[pidAT+idx]),
                             Real3D(vx[pidAT+idx], vy[pidAT+idx], vz[pidAT+idx]),
                             Real3D(fx[pidAT+idx], fy[pidAT+idx], fz[pidAT+idx]), 1, 1, 1])
        
    tuples.append(tmptuple) # add tuples
    idx = idx+4

# add particles
system.storage.addParticles(allParticles, "id", "pos", "v", "f", "type", "mass", "adrat")

# add tuples
ftpl = espresso.FixedTupleList(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuples(ftpl)

system.storage.decompose()

# AdResS Verlet List
vl = espresso.VerletListAdress(system, cutoff=rc, dEx=12, dHy=2.5, adrCenter=[18.42225, 18.42225, 18.42225])

# Non-bonded potentials
# Tabulated WCA between AT and tabulated Morse between CG particles
interNB = espresso.interaction.VerletListAdressTabulated(vl, ftpl)  
potWCA = espresso.interaction.Tabulated(itype=2, filename=tabWCA, cutoff=rca)
potMorse = espresso.interaction.Tabulated(itype=2, filename=tabMorse, cutoff=rc)
interNB.setPotential(type1=0, type2=0, potential=potMorse)
interNB.setPotential(type1=1, type2=1, potential=potWCA)
system.addInteraction(interNB)


# FENE bonds between AT particles
fpl = espresso.FixedPairListAdress(system.storage)
#fpl.addBonds(bonds)
potFENE = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espresso.interaction.FixedPairListFENE(system, fpl, potFENE)
system.addInteraction(interFENE)


# setup integrator
integrator = espresso.integrator.VelocityVerletAdress(system)
integrator.dt = timestep

if(nvt):
  langevin = espresso.integrator.Langevin(system)
  langevin.gamma = 0.5
  langevin.temperature = 1.0
  integrator.langevin = langevin
  integrator.dt = timestep

print ''
print 'number of AT particles =', num_particles
print 'number of CG particles =', num_particlesCG
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'nvt =', nvt
print 'steps =', steps
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)
pressureTensor = espresso.analysis.PressureTensor(system)

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
for i in range(1, intervals + 1):
  integrator.run(nsteps)
  step = nsteps * i
  T = temperature.compute()
  P = pressure.compute()
  Pij = pressureTensor.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interNB.computeEnergy()
  Eb = interFENE.computeEnergy()
  sys.stdout.write(fmt % (step, T, P, Pij[3], Ek + Ep + Eb, Ep, Eb, Ek))
end_time = time.clock()


timers.show(integrator.getTimers(), precision=3)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))

