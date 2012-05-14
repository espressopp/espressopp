#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import time
import espresso
import MPI
import logging
from espresso import Real3D, Int3D
from espresso.tools.convert import gromacs
from espresso.tools import decomp
from espresso.tools import timers


# integration steps, timestep, intervals
steps = 10000
timestep = 0.001
intervals = 10*2

# cutoffs, skin
rc = 2.31 # CG cutoff, Morse
rca = 1.122462048309373 # AT cutoff (2^(1/6)), WCA
skin = 0.4

# AdResS setup
ex_size = 12.0
hy_size = 2.5
adrCenter = [18.42225, 18.42225, 18.42225]

# Langevin thermostat setup
gamma = 0.5
temp = 1.0

# PDB output
pdb_dir = "./pdb/"
writepdb = False

# GROMACS configuration file
init_conf = "conf.gro"



#### helper functions ####

# makes a pdb from a list of particles
def make_pdb(system, num_particles, size, unfold=False, outfile="output", step=0):
    step = ("%04d" % step)  
    pdb = open(pdb_dir+outfile+step+".pdb", 'w')
    mol = 0
    moltmp = 0
    for i in range(num_particles):
        part = system.storage.getParticle(i)
        moltmp = moltmp+1
        pos = part.pos
        
        if (unfold):
            # unfold position
            for dir in range(3):
                pos[dir] += part.imageBox[dir] * size[dir]
        
        #ATOM      1  N   ALA A   2      12.953  66.007  46.200  1.00 39.81           N  
        pdb.write("%-6s%5d  %-3s %3s %1s%4s %11.3f %7.3f %7.3f  0.00 00.00 %11s  \n" %
                  ("ATOM", part.id, "FE", "UNX", "F", mol, pos[0], pos[1], pos[2], "T00"+str(part.type)))
        if moltmp == 4:
            mol = mol+1
            moltmp = 0
    pdb.close()

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

##########################


# create tabulated morse potential used for CG interactions
tabMorse = "pot-morse.txt"
potMorse = espresso.interaction.Morse(epsilon=0.105, alpha=2.4, rMin=rc, cutoff=rc, shift="auto")
writeTabFile(potMorse, tabMorse, N=512, low=0.005, high=rc)


# read GROMACS file with inital configuration
x, y, z, vx, vy, vz, Lx, Ly, Lz = gromacs.read(init_conf)
num_particles = len(x)

sys.stdout.write('Setting up simulation ...\n')

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

# AdResS domain decomposition
system.storage = espresso.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)


# prepare CG and AT particles with properties
props = ['id', 'pos', 'v', 'type', 'mass', 'adrat']
particles, tuples, bonds = [], [], []
for i in range(0, num_particles, 5):
    
    # CG particle
    particles.append([i, Real3D(x[i], y[i], z[i]),
                      Real3D(vx[i], vy[i], vz[i]), 0, 4.0, 0])
    
    # AT particles
    for j in range(1, 5):
        p = i+j
        
        particles.append([p, Real3D(x[p], y[p], z[p]),
                      Real3D(vx[p], vy[p], vz[p]), 1, 1.0, 1])
        
        # prepare bonds
        for k in range(j+1, 5): 
            bonds.append([i+j, i+k])
        
    tuples.append([i, i+1, i+2, i+3, i+4])

# add particles
sys.stdout.write('Adding particles ...\n')
system.storage.addParticles(particles, *props)

# add tuples
ftpl = espresso.FixedTupleList(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuples(ftpl)
        
# add bonds between AT particles
fpl = espresso.FixedPairListAdress(system.storage, ftpl)
fpl.addBonds(bonds)

# decompose after adding tuples and bonds
sys.stdout.write('Added tuples and bonds, decomposing now ...\n')
system.storage.decompose()


# AdResS Verlet list
vl = espresso.VerletListAdress(system, cutoff=rc+skin, adrcut=rc+skin,
                                dEx=ex_size, dHy=hy_size,
                                adrCenter=adrCenter)

# non-bonded potentials
# LJ Capped WCA between AT and tabulated Morse between CG particles
# exclusions not needed, since there are no intramolecular interactions 
interNB = espresso.interaction.VerletListAdressLennardJonesCapped(vl, ftpl)
potWCA  = espresso.interaction.LennardJonesCapped(epsilon=1.0, sigma=1.0, shift=True, caprad=0.27, cutoff=rca)
potMorse = espresso.interaction.Tabulated(itype=2, filename=tabMorse, cutoff=rc) # CG
interNB.setPotentialAT(type1=1, type2=1, potential=potWCA) # AT
interNB.setPotentialCG(type1=0, type2=0, potential=potMorse) # CG
system.addInteraction(interNB)

# bonded potentials
# FENE and LJ potential between AT particles
potFENE = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
potLJ = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift=True, cutoff=rca)
interFENE = espresso.interaction.FixedPairListFENE(system, fpl, potFENE)
interLJ = espresso.interaction.FixedPairListLennardJones(system, fpl, potLJ)
system.addInteraction(interFENE)
system.addInteraction(interLJ)


# AdResS integrator
integrator = espresso.integrator.VelocityVerletAdress(system)
integrator.dt = timestep

# Langevin thermostat
langevin = espresso.integrator.Langevin(system)
langevin.gamma = gamma
langevin.temperature = temp
integrator.langevin = langevin
integrator.dt = timestep


print ''
print 'steps =', steps
print 'dt =', integrator.dt
print 'number of CG+AT particles =', num_particles
print 'rc CG =', rc
print 'rc AT =', rca
print 'skin =', system.skin
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
temperature = espresso.analysis.Temperature(system)

fmt = '%5d %8.4f %12.3f %12.3f %12.3f %12.3f\n'

T = temperature.compute()
Ek = 0.5 * T * (3 * num_particles)
Ep = interNB.computeEnergy()
Eb = interFENE.computeEnergy()
sys.stdout.write(' step        T       etotal   epotential      ebonded     ekinetic\n')
sys.stdout.write(fmt % (0, T, Ek + Ep + Eb, Ep, Eb, Ek))

start_time = time.clock()
nsteps = steps / intervals
for s in range(1, intervals + 1):
  integrator.run(nsteps)
  step = nsteps * s
  T = temperature.compute()
  Ek = 0.5 * T * (3 * num_particles)
  Ep = interNB.computeEnergy()
  Eb = interFENE.computeEnergy()
  sys.stdout.write(fmt % (step, T, Ek + Ep + Eb, Ep, Eb, Ek))
  system.storage.decompose()

  if (writepdb):
      make_pdb(system, num_particles, size, unfold=False, outfile="atomsADR", step=s) # AT particles

end_time = time.clock()

timers.show(integrator.getTimers(), precision=3)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))

