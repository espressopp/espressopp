#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################################################
#  ESPResSo++                                                             #
#  Test script for Tabulated potentials (LJ, FENE, Cosine)                #
#                                                                         #
###########################################################################

import sys
import time
import espresso
import MPI
import math
import logging
import os
from espresso import Real3D, Int3D
from espresso.tools.convert import lammps, gromacs
from espresso.tools import decomp
from espresso.tools import timers

# simulation parameters (nvt = False is nve)
steps = 100
rc = 1.12
skin = 0.3                                 # skin for Verlet lists
nvt = True
timestep = 0.01
spline  = 2                                # spline interpolation type (1, 2, 3)

conffile = 'polymer_melt.start' # file with inital configuration
tabfileLJ = "pot-lj.txt"
tabfileFENE = "pot-fene.txt"
tabfileCosine = "pot-cosine.txt"


######################################################################
##  IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE  ##
######################################################################

bonds, angles, x, y, z, Lx, Ly, Lz = lammps.read(conffile)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

print '\n-- Tabulated Potentials Test --\n'
print 'Steps: %3s' % steps
print 'Particles: %3s' % num_particles
print 'Cutoff: %3s' % rc
print 'Density = %.4f' % (density)
print 'dt =', timestep
print 'Skin =', skin
print 'nvt =', nvt

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

# compute the number of cells on each node
def calcNumberCells(size, nodes, cutoff):
    ncells = 1
    while size / (ncells * nodes) >= cutoff:
       ncells = ncells + 1
    return ncells - 1


# write the tabulated potential files
potLJ  = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift=False, cutoff=rc)
potFENE = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
potCosine = espresso.interaction.Cosine(K=1.5, theta0=3.1415926)

print 'Generating potential files ... (%2s, %2s, %2s)\n' % (tabfileLJ, tabfileFENE, tabfileCosine)
writeTabFile(potLJ, tabfileLJ, N=257, low=0.01, high=potLJ.cutoff)
writeTabFile(potFENE, tabfileFENE, N=257, low=0.0001, high=1.49)
writeTabFile(potCosine, tabfileCosine, N=257, low=0.0001, high=3.14, body=3)

potTabLJ = espresso.interaction.Tabulated(itype=spline, filename=tabfileLJ, cutoff=rc)
potTabFENE = espresso.interaction.Tabulated(itype=spline, filename=tabfileFENE)
potTabCosine = espresso.interaction.TabulatedAngular(itype=spline, filename = tabfileCosine)

# repeat simulation twice, with and without tabulated potential
for tabulation in [True, False]:
    if tabulation: w = ''
    else: w = 'out'
    print 'Running simulation with%0s tabulated potentials' % w
        
    print 'Setting up ...'
    system = espresso.System()
    system.rng = espresso.esutil.RNG(54321)
    system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
    system.skin = skin
        
    comm = MPI.COMM_WORLD
        
    nodeGrid = decomp.nodeGrid(comm.size)
    cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
    #nodeGrid = Int3D(1, 1, comm.size)
    #cellGrid = Int3D(
        #calcNumberCells(size[0], nodeGrid[0], rc),
        #calcNumberCells(size[1], nodeGrid[1], rc),
        #calcNumberCells(size[2], nodeGrid[2], rc)
        #)
    system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
        
    # add particles to the system and then decompose
    for pid in range(num_particles):
        system.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
    system.storage.decompose()
        
        
    # Lennard-Jones with Verlet list
    vl = espresso.VerletList(system, cutoff = rc + system.skin)
    if tabulation:
        interLJ = espresso.interaction.VerletListTabulated(vl)
        interLJ.setPotential(type1=0, type2=0, potential=potTabLJ)
    else:
        interLJ = espresso.interaction.VerletListLennardJones(vl)
        interLJ.setPotential(type1=0, type2=0, potential=potLJ)
    system.addInteraction(interLJ)
        
        
    # FENE bonds with Fixed Pair List
    fpl = espresso.FixedPairList(system.storage)
    fpl.addBonds(bonds)
    if tabulation:
        interFENE = espresso.interaction.FixedPairListTabulated(system, fpl, potTabFENE)
        #interFENE.setPotential(type1=0, type2=0, potential=potTabFENE) # no longer needed
    else:
        interFENE = espresso.interaction.FixedPairListFENE(system, fpl, potFENE)
        #interFENE.setPotential(type1=0, type2=0, potential=potFENE)
    system.addInteraction(interFENE)
        
        
    # Cosine with Fixed Triple List
    ftl = espresso.FixedTripleList(system.storage)
    ftl.addTriples(angles)
    if tabulation:
        interCosine = espresso.interaction.FixedTripleListTabulatedAngular(system, ftl, potTabCosine)
        #interCosine.setPotential(type1=0, type2=0, potential=potTabCosine)
    else:
        interCosine = espresso.interaction.FixedTripleListCosine(system, ftl, potCosine)
        #interCosine.setPotential(type1=0, type2=0, potential=potCosine)
    system.addInteraction(interCosine)
        
        
        
    # integrator
    integrator = espresso.integrator.VelocityVerlet(system)
    integrator.dt = timestep
        
    if(nvt):
        langevin = espresso.integrator.LangevinThermostat(system)
        langevin.gamma = 1.0
        langevin.temperature = 1.0
        integrator.addExtension(langevin)
        
        
    # analysis
    configurations = espresso.analysis.Configurations(system)
    configurations.gather()
    temperature = espresso.analysis.Temperature(system)
    pressure = espresso.analysis.Pressure(system)
    pressureTensor = espresso.analysis.PressureTensor(system)
        
    fmt = '%5d %8.4f %10.5f %8.5f %12.3f %12.3f %12.3f %12.3f %12.3f\n'
        
    T = temperature.compute()
    P = pressure.compute()
    Pij = pressureTensor.compute()
    Ek = 0.5 * T * (3 * num_particles)
    Ep = interLJ.computeEnergy()
    Eb = interFENE.computeEnergy()
    Ea = interCosine.computeEnergy()
    Etotal = Ek + Ep + Eb + Ea
    sys.stdout.write(' step     T          P       Pxy        etotal      ekinetic      epair        ebond       eangle\n')
    sys.stdout.write(fmt % (0, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))
        
    start_time = time.clock()
    integrator.run(steps)
        
    T = temperature.compute()
    P = pressure.compute()
    Pij = pressureTensor.compute()
    Ek = 0.5 * T * (3 * num_particles)
    Ep = interLJ.computeEnergy()
    Eb = interFENE.computeEnergy()
    Ea = interCosine.computeEnergy()
    Etotal = Ek + Ep + Eb + Ea
    sys.stdout.write(fmt % (steps, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea))
    sys.stdout.write('\n')
        
    end_time = time.clock()
    timers.show(integrator.getTimers(), precision=2)

#os.system('rm '+tabfileLJ+' '+tabfileFENE+' '+tabfileCosine)
print '\nDone.'
