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

###########################################################################
#  ESPResSo++ Python script for H-AdResS Water                            #
#  simulation  based on Gromacs topology                                  #
###########################################################################
	
import math
import sys
import time
import espressopp
import mpi4py.MPI as MPI
import logging
from espressopp import Real3D, Int3D
from espressopp.tools import gromacs
from espressopp.tools import decomp
from espressopp.tools import timers

# This example reads in a gromacs water system (SPC/Fw) treated with reaction field. See the corresponding gromacs grompp.mdp paramter file.
# Output of gromacs energies and esp energies should be the same
# For H-Adress, special interactions and domain decomposition have to be defined. The gromacs parser has an option to create Adress interactions instead of standard ones
# In the current implementation only one type of atomistic potential can be set for each interaction(template). This makes it necessarry to create two interaction templates (one for coulomb, one for lennard-jones) and leave the coarse-grained interaction unset in one of them.

# simulation parameters (nvt = False is nve)
steps = 100000
check = steps/1
timestep = 0.0005

# parameters to convert GROMACS tabulated potential file
sigma = 1.0
epsilon = 1.0
c6 = 1.0
c12 = 1.0

# H-AdResS
rc = 1.3  # cutoff coarse-grained potential
rca = 0.9  # cutoff atomistic potential
skin = 0.2

# parameters for size of AdResS dimensions
ex_size = 1.5
hy_size = 2.0

# GROMACS setup files
grofile = "conf.gro"
topfile = "topol.top"


# this calls the gromacs parser for processing the top file (and included files) and the conf file
# The variables at the beginning defaults, types, etc... can be found by calling
# gromacs.read(grofile,topfile) without return values. It then prints out the variables to be unpacked
defaults, types, atomtypes, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, exclusions, x, y, z, vx, vy, vz, resname, resid, Lx, Ly, Lz =gromacs.read(grofile,topfile)

# this is an equilibrated configuration!
dummy1, dummy2, x, y, z, vx, vy, vz, dummy3, dummy4, dummy5 = espressopp.tools.readxyz("equilibrated_conf.xyz")

# particles, geometry, density
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

# set up the system
sys.stdout.write('Setting up simulation ...\n')
system = espressopp.System()

# random number generator
xs = time.time()
seed = int(xs % int(xs) * 10000000000)
rng = espressopp.esutil.RNG()
rng.seed(seed)
system.rng = rng

# boundary conditions
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

# communication, storage and cell/node grid
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size,size,rc,skin)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

# create a verlet list
verletlist = espressopp.VerletListAdress(system, cutoff=rc, adrcut=rc,
                                dEx=ex_size, dHy=hy_size,
                                adrCenter=[Lx/2, Ly/2, Lz/2])

# add particles to the system and then decompose
props = ['id', 'pos', 'v', 'f', 'type', 'mass', 'q', 'adrat']
allParticlesAT = []
allParticles = []
tuples = []
# prepare AT particles
for pid in range(num_particles):
    part = [pid + 1, Real3D(x[pid], y[pid], z[pid]),
            Real3D(vx[pid],vy[pid], vz[pid]), Real3D(0, 0, 0),
            types[pid], masses[pid], charges[pid], 1]
    allParticlesAT.append(part)

num_particlesCG = len(x)/3
typeCG=0
# create CG particles
for pidCG in range(num_particlesCG):
    # we put CG molecule in first atom, later CG molecules will be positioned in the center
    #cmp = espressopp.tools.AdressSetCG(3, pidCG, allParticlesAT)

    # Preparation of tuples (tuples define, which atoms belong to which CG molecules)
    tmptuple = [pidCG+num_particles+1]
    for pidAT2 in range(3):
        pid = pidCG*3+pidAT2
        tmptuple.append((allParticlesAT[pid])[0])
    firsParticleId=tmptuple[1]
    cmp=allParticlesAT[firsParticleId-1][1]

    typeCG=max(types)+1
    # append CG particles
    allParticles.append([pidCG+num_particles+1, # CG particle has to be added first!
                         Real3D(cmp[0], cmp[1], cmp[2]), # pos
                         Real3D(0, 0, 0), # vel
                         Real3D(0, 0, 0), # force
                         typeCG, 18.0154, 0.0, 0]) # type, mass, q, is not AT particle
    # append AT particles
    for pidAT in range(3):
        pid = pidCG*3+pidAT
        allParticles.append([(allParticlesAT[pid])[0], # now the AT particles can be added
                            (allParticlesAT[pid])[1], # pos
                            (allParticlesAT[pid])[2], # vel
                            (allParticlesAT[pid])[3], # force
                            (allParticlesAT[pid])[4], # type
                            (allParticlesAT[pid])[5], # mass
                            (allParticlesAT[pid])[6], # q
                            (allParticlesAT[pid])[7]]) # is AT particle
    # append tuple to tuplelist
    tuples.append(tmptuple)

system.storage.addParticles(allParticles, *props)


# create FixedTupleList object and add the tuples
ftpl = espressopp.FixedTupleListAdress(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)
system.storage.decompose()

# set up LJ interaction according to the parameters read from the .top file
ljinteraction=gromacs.setLennardJonesInteractions(system, defaults, atomtypeparameters, verletlist,rca, hadress=True, ftpl=ftpl)

# set up angle interactions according to the parameters read from the .top file
angleinteractions=gromacs.setAngleInteractionsAdress(system, angletypes, angletypeparams, ftpl)

# set up coulomb interactions according to the parameters read from the .top file
# !! Warning: this only works for reaction-field now!
qq_interactions=gromacs.setCoulombInteractions(system, verletlist, rca, types, epsilon1=1, epsilon2=80, kappa=0, hadress=True, ftpl=ftpl)

# load CG interaction from table
fe="table_CG_CG.tab"
gromacs.convertTable("table_CG_CG.xvg", fe, 1, 1, 1, 1)
potCG = espressopp.interaction.Tabulated(itype=3, filename=fe, cutoff=rca) # CG

# set the CG potential. There are two non-bonded interactions, we pick only the first one
for n in range(system.getNumberOfInteractions()):
    interaction=system.getInteraction(n)
    if interaction.bondType() == espressopp.interaction.Nonbonded:
	print "Setting CG interaction", typeCG
	interaction.setPotentialCG(type1=typeCG, type2=typeCG, potential=potCG)
	break

# set up bonded interactions according to the parameters read from the .top file
bondedinteractions=gromacs.setBondedInteractionsAdress(system, bondtypes, bondtypeparams, ftpl)

# exlusions, i.e. pairs of atoms not considered for the non-bonded part. Those are defined either by bonds which automatically generate an exclusion. Or by the nregxcl variable
verletlist.exclude(exclusions)

# add VelocityVerlet Integrator
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.dt = timestep

# add Langevin Thermostat
langevin = espressopp.integrator.LangevinThermostat(system)
langevin.gamma = 2.0
langevin.temperature = 2.4942 # kT in gromacs units
langevin.adress = True
integrator.addExtension(langevin)

# add AdResS
adress = espressopp.integrator.Adress(system,verletlist,ftpl)
integrator.addExtension(adress)

# distribute atoms and CG molecules according to AdResS domain decomposition, place CG molecules in the center of mass
espressopp.tools.AdressDecomp(system, integrator)

# print simulation parameters
print ''
print 'number of particles =', num_particles
print 'density = %.4f' % (density)
print 'rc =', rc
print 'dt =', integrator.dt
print 'skin =', system.skin
print 'steps =', steps
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

# analysis
configurations = espressopp.analysis.Configurations(system)
configurations.gather()
temperature = espressopp.analysis.Temperature(system)
pressure = espressopp.analysis.Pressure(system)
pressureTensor = espressopp.analysis.PressureTensor(system)

print "i*timestep, T, Eb, EAng, ELj, EQQ, Ek, Etotal"
fmt='%5.5f %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8f\n'

start_time = time.clock()
outfile = open("esp.dat", "w")

# write a snapshot of the system
espressopp.tools.psfwrite("system.psf", system, typenames={0:'H', 1:'O', 2:'CG'})
espressopp.tools.pdbwrite("system.pdb", system, append=False, typenames={0:'H', 1:'O', 2:'CG'})

for i in range(check):

    T = temperature.compute()
    P = pressure.compute()
    Eb = 0
    EAng = 0
    for bd in bondedinteractions.values(): Eb+=bd.computeEnergy()
    for ang in angleinteractions.values(): EAng+=ang.computeEnergy()
    ELj= ljinteraction.computeEnergy()
    EQQ= qq_interactions.computeEnergy()
    Ek = 0.5 * T * (3 * num_particles)
    Etotal = Ek+Eb+EAng+EQQ+ELj
    outfile.write(fmt%(i*steps/check*timestep, T, Eb, EAng, ELj, EQQ, Ek, Etotal))
    print (fmt%(i*steps/check*timestep, T, Eb, EAng, ELj, EQQ, Ek, Etotal))

    integrator.run(steps/check) # print out every steps/check steps


# simulation information
end_time = time.clock()
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))
