#!/usr/bin/env python
# -*- coding: utf-8 -*-



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
steps = 10000
rc = 2.31 # CG cutoff, Morse
rca = 1.122462048309373 # AT cutoff (2^(1/6)), WCA
skin = 0.8
nvt = True
timestep = 0.001
intervals = 10

ex_size = 12.0
hy_size = 2.5
pdb_dir = "./pdbs-test/"
writepdb = False

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
                pos[dir] += part.image[dir] * size[dir]
        
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
            force /= r
        else: # this is for 3- and 4-body potentials
            force = pot.computeForce(r)
        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))
     
    outfile.close()


# capped LJ
tabWCA = "pot-wca.txt"
potLJ  = espresso.interaction.LennardJonesCapped(epsilon=1.0, sigma=1.0, shift=True, caprad=0.27, cutoff=rca)
writeTabFile(potLJ, tabWCA, N=512, low=0.005, high=rca)

#tabFENE = "pot-fene.txt"
#potFENE  = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
#writeTabFile(potFENE, tabFENE, N=512, low=0.005, high=2.0)

tabMorse = "pot-morse.txt"
potMorse = espresso.interaction.Morse(epsilon=0.105, alpha=2.4, rMin=rc, cutoff=rc, shift="auto")
writeTabFile(potMorse, tabMorse, N=512, low=0.005, high=rc)



# read ESPResSo configuration file 
Lx, Ly, Lz, x, y, z, type, q, vx, vy, vz, fx, fy, fz, bonds = espresso_old.read("ex_j0.txt")
num_particlesCG = 5001 # number of VP/CG particles
#num_particles = len(x) - num_particlesCG  # 20004 = 25005 - 5001 
num_particles = len(x) # 20004


#Lx, Ly, Lz = 45, 45, 45

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

# AdResS domain decomposition
system.storage = espresso.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)


## add CG and AT particles
#allParticles = []
#tuples = []
#idx = 0 # just an index
#for pidCG in range(num_particles, num_particles+num_particlesCG):
#    #print pidCG,
#    tmptuple = [pidCG] # first pid in tuple is pid of CG/VP particle
#    allParticles.append([pidCG,
#                         Real3D(x[pidCG], y[pidCG], z[pidCG]),
#                         Real3D(vx[pidCG], vy[pidCG], vz[pidCG]),
#                         Real3D(fx[pidCG], fy[pidCG], fz[pidCG]), 0, 4, 0])
#    
#    for pidAT in range(4): # each CG has 4 AT particles
#        #print pidAT+idx,
#        tmptuple.append(pidAT+idx) # other pids in tuple are pids of AT particles
#        allParticles.append([pidAT+idx, Real3D(x[pidAT+idx], y[pidAT+idx], z[pidAT+idx]),
#                             Real3D(vx[pidAT+idx], vy[pidAT+idx], vz[pidAT+idx]),
#                             Real3D(fx[pidAT+idx], fy[pidAT+idx], fz[pidAT+idx]), 1, 1, 1])
#        
#    tuples.append(tmptuple) # add tuples
#    idx += 4



# add AT particles
allParticlesAT = []
allParticles = []
tuples = []
for pidAT in range(num_particles):
    #print pidAT,
    allParticlesAT.append([pidAT, # add these particles just temporarly! 
                         Real3D(x[pidAT], y[pidAT], z[pidAT]),
                         Real3D(vx[pidAT], vy[pidAT], vz[pidAT]),
                         Real3D(fx[pidAT], fy[pidAT], fz[pidAT]), 1, 1, 1])

# create CG particles from center of mass
for pidCG in range(num_particlesCG):
    cmp = [0,0,0]
    cmv = [0,0,0]
    tmptuple = [pidCG+num_particles]
    for pidAT in range(4):
        tmptuple.append(pidCG*4+pidAT)
        pos = (allParticlesAT[pidCG*4+pidAT])[1]
        vel = (allParticlesAT[pidCG*4+pidAT])[2]
        cmp[0] += pos[0]
        cmp[1] += pos[1]
        cmp[2] += pos[2]
        cmv[0] += vel[0]
        cmv[1] += vel[1]
        cmv[2] += vel[2]
    cmp[0] /= 4
    cmp[1] /= 4
    cmp[2] /= 4
    cmv[0] /= 4
    cmv[1] /= 4
    cmv[2] /= 4
    allParticles.append([pidCG+num_particles, # CG particle has to bo added first!
                         Real3D(cmp[0], cmp[1], cmp[2]),
                         Real3D(cmv[0], cmv[1], cmv[2]),
                         Real3D(0, 0, 0), 0, 4, 0])
    for pidAT in range(4): 
        allParticles.append([pidCG*4+pidAT, # now the AT particles can be added
                         Real3D(x[pidCG*4+pidAT], y[pidCG*4+pidAT], z[pidCG*4+pidAT]),
                         Real3D(vx[pidCG*4+pidAT], vy[pidCG*4+pidAT], vz[pidCG*4+pidAT]),
                         Real3D(fx[pidCG*4+pidAT], fy[pidCG*4+pidAT], fz[pidCG*4+pidAT]), 1, 1, 1])
        
    tuples.append(tmptuple)
    
#print allParticles
#print tuples
#exit()

# add particles
system.storage.addParticles(allParticles, "id", "pos", "v", "f", "type", "mass", "adrat")

#exit()

# add tuples
ftpl = espresso.FixedTupleList(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuples(ftpl)

# add bonds between AT particles
fpl = espresso.FixedPairListAdress(system.storage, ftpl)
fpl.addBonds(bonds)

# decompose after adding tuples and bonds
print "\n\nAdded tuples and bonds, decomposing now ..." 
system.storage.decompose()


# AdResS Verlet list
vl = espresso.VerletListAdress(system, cutoff=rc+skin, dEx=ex_size, dHy=hy_size, adrCenter=[18.42225, 18.42225, 18.42225])

# non-bonded potentials
# Tabulated WCA between AT and tabulated Morse between CG particles
interNB = espresso.interaction.VerletListAdressTabulated(vl, ftpl)
potWCA = espresso.interaction.Tabulated(itype=2, filename=tabWCA, cutoff=rca) # AT
potMorse = espresso.interaction.Tabulated(itype=2, filename=tabMorse, cutoff=rc) # CG
interNB.setPotential(type1=0, type2=0, potential=potMorse) # CG
interNB.setPotential(type1=1, type2=1, potential=potWCA) # AT
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

  if (writepdb):
      #  for i in range(num_particles-1, num_particles+num_particlesCG): # CG particles
      #num_particles = num_particles+num_particlesCG # all particles
      make_pdb(system, num_particles, size, unfold=False, outfile="atomsADR", step=s) # AT particles

end_time = time.clock()

#start_time = time.clock()
#integrator.run(steps)
#T = temperature.compute()
#P = pressure.compute()
#Pij = pressureTensor.compute()
#Ek = 0.5 * T * (3 * num_particles)
#Ep = interNB.computeEnergy()
#Eb = interFENE.computeEnergy()
#end_time = time.clock()
#sys.stdout.write(fmt % (steps, T, P, Pij[3], Ek + Ep + Eb, Ep, Eb, Ek))

timers.show(integrator.getTimers(), precision=3)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))

