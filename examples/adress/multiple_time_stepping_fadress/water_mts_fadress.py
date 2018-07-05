#!/usr/bin/env python
#  Copyright (C) 2018
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

############################################################################################
#                                                                                                                                                                   #
#  ESPResSo++ Python script for an F-AdResS simulation of water/IBI potential using Multiple Time Stepping  #
#                                                                                                                                                                    #
#############################################################################################

import time
import espressopp
import mpi4py.MPI as MPI
from espressopp import Real3D
from espressopp.tools.convert import gromacs
from espressopp.tools import decomp

# Performs an F-AdResS adaptive resolution simulation using a RESPA multiple time stepping scheme.
# All atomistic interactions (both bonded and non-bonded) are integrated on the fast time scale, while coarse-grained
# interactions are calculated on the slow time scale. The atomistic force field is SPC/Fw with reaction field
# while the coarse-grained interaction potential is a potential derived from iterative Boltzmann inversion.
print 'Performing an ESPResSo++ F-AdResS simulation with multiple time stepping.\n'

########################################################
#  1. specification of the system setup and simulation parameters  #
########################################################

print 'Specifying simulation parameters...'
totaltime = 0.25 # total runtime (in ps)
printfreq = 0.005 # frequency (in ps) to print outputs
timestep = 0.0005 # short timestep
longstep = 0.0025 # long timestep

multistep = int(longstep/timestep) # short steps per single long step
steps = int(totaltime/longstep) # total number of long steps
checkfreq = int(printfreq/longstep) # number of long time steps between printed outputs/checks
check = steps/checkfreq # number of intervals to run the integrator between checks

interaction_cutoff_cg = 1.0  # cutoff coarse-grained potential
interaction_cutoff_at = 1.0  # cutoff atomistic potential
interaction_cutoff_verletlist_at = 1.25 # cutoff atomistic potential in verlet list cutoff (larger than normal atomistic interaction cutoff, since verletlist is built using CG particles)
skin = 0.05 # verlet list skin

ex_size = 2.2 # radius of atomistic region
hy_size = 1.0 # width of hybrid region

correction_on_slow_timescale = False # apply corrections on slow timescale?

######################
#  2. read in coordinates  #
######################

print 'Reading in coordinates... (not equilibrated!)'
# call the gromacs parser for processing the top file (and included files) and the conf file (note that this configuration is most likely not equilibrated)
defaults, types, atomtypes, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, exclusions, x, y, z, vx, vy, vz, resname, resid, Lx, Ly, Lz = gromacs.read("conf.gro", "topol.top")

# particles, geometry, density
num_particles_AT = len(x)
num_particles_CG = len(x)/3
density = num_particles_AT / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

#####################
#  3. set up the system  #
#####################

print 'Setting up system...'
# system, boundary conditions, skin, communicator, node & cell grids
system = espressopp.System()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, interaction_cutoff_verletlist_at, skin)

# random number generator
rng = espressopp.esutil.RNG()
rng.seed(42)
system.rng = rng

# AdResS domain decomposition
system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

##########################
#  4. add particles to system  #
##########################

print 'Adding particles and tuples...'
props = ['id', 'pos', 'v', 'f', 'type', 'mass', 'q', 'adrat']
allParticlesAT = []
allParticles = []
tuples = []

# prepare atomistic particles (add these particles here just temporarily)
for pid in range(num_particles_AT):
    part = [pid + 1, Real3D(x[pid], y[pid], z[pid]),
            Real3D(vx[pid],vy[pid], vz[pid]), Real3D(0, 0, 0),
            types[pid], masses[pid], charges[pid], 1]
    allParticlesAT.append(part)

# create coarse-grained particles
for pidCG in range(num_particles_CG):

    # we put CG molecule in first atom, later CG molecules will be positioned in the center
    # note that we put the CG molecule at the first atom's position. Later the CG molecule will be positioned in the center of mass of all it's atoms
    tmptuple = [pidCG+num_particles_AT+1]
    for pidAT in range(3):
        pid = pidCG*3+pidAT
        tmptuple.append((allParticlesAT[pid])[0])
    firsParticleId=tmptuple[1]
    cmp=allParticlesAT[firsParticleId-1][1]
    typeCG=max(types)+1

    # append coarse-grained particles
    allParticles.append([pidCG+num_particles_AT+1,
                         Real3D(cmp[0], cmp[1], cmp[2]), # pos
                         Real3D(0, 0, 0), # vel
                         Real3D(0, 0, 0), # force
                         typeCG, 18.0154, 0.0, 0]) # type, mass, q, is not AT particle
    # append atomistic particles
    for pidAT in range(3):
        pid = pidCG*3+pidAT
        allParticles.append([(allParticlesAT[pid])[0],
                            (allParticlesAT[pid])[1], # pos
                            (allParticlesAT[pid])[2], # vel
                            (allParticlesAT[pid])[3], # force
                            (allParticlesAT[pid])[4], # type
                            (allParticlesAT[pid])[5], # mass
                            (allParticlesAT[pid])[6], # q
                            (allParticlesAT[pid])[7]]) # is AT particle
    # append tuple to tuplelist
    tuples.append(tmptuple)

# add particles to system
system.storage.addParticles(allParticles, *props)

# create FixedTupleList object and add the tuples
ftpl = espressopp.FixedTupleListAdress(system.storage)
ftpl.addTuples(tuples)
system.storage.setFixedTuplesAdress(ftpl)

# decompose to distribute the particles correctly before adding interactions
system.storage.decompose()

###########################################
#  5. set up structure, interactions, and force field  #
###########################################

print 'Setting up interactions and force field...'
# create the adaptive resolution verlet list
verletlist = espressopp.VerletListAdress(system, cutoff=interaction_cutoff_cg, adrcut=interaction_cutoff_verletlist_at,
                                         dEx=ex_size, dHy=hy_size, adrCenter=[Lx/2, Ly/2, Lz/2])

# set up atomistic nonbonded interactions. This interaction template handles both the Lannard Jones and the Reaction Field within AdResS
interNBat = espressopp.interaction.VerletListAdressATLenJonesReacFieldGen(verletlist, ftpl)
potLJ  = espressopp.interaction.LennardJones(epsilon=0.650299305951, sigma=0.316549165245, shift='auto', cutoff=interaction_cutoff_at)
interNBat.setPotential1(type1=1, type2=1, potential=potLJ)
potQQ  = espressopp.interaction.ReactionFieldGeneralized(prefactor=138.935485, kappa=0.0, epsilon1=1.0, epsilon2=80.0, cutoff=interaction_cutoff_at, shift="auto")
interNBat.setPotential2(type1=1, type2=1, potential=potQQ)
interNBat.setPotential2(type1=1, type2=0, potential=potQQ)
interNBat.setPotential2(type1=0, type2=0, potential=potQQ)
system.addInteraction(interNBat)

# set up coarse-grained potential (an iterative Boltzmann Inversion potential)
potCG = espressopp.interaction.Tabulated(itype=3, filename="table_ibi.dat", cutoff=interaction_cutoff_cg)
interNBcg = espressopp.interaction.VerletListAdressCGTabulated(verletlist, ftpl)
interNBcg.setPotential(type1=typeCG, type2=typeCG, potential=potCG)
system.addInteraction(interNBcg)

# set up angle interactions
angleinteractions=gromacs.setAngleInteractionsAdress(system, angletypes, angletypeparams, ftpl)

# set up bonded interactions
bondedinteractions=gromacs.setBondedInteractionsAdress(system, bondtypes, bondtypeparams, ftpl)

# verletlist exlusions, i.e. pairs of atoms not considered for the non-bonded force calculation
verletlist.exclude(exclusions)

#########################################
#  6. set up integration scheme and corrections  #
#########################################

print 'Setting up integration scheme (RESPA velocity verlet integrator)...'
# set up the RESPA VelocityVerlet Integrator
integrator = espressopp.integrator.VelocityVerletRESPA(system)
integrator.dt = timestep
integrator.multistep = multistep

# add AdResS
adress = espressopp.integrator.Adress(system, verletlist, ftpl, multistep=multistep)
integrator.addExtension(adress)

# add Langevin Thermostat
langevin = espressopp.integrator.LangevinThermostat(system)
langevin.gamma = 10.0
langevin.temperature = 2.494338 # kT in gromacs units
langevin.adress = True
integrator.addExtension(langevin)

# add Thermodynamic Force
thdforce = espressopp.integrator.TDforce(system, verletlist, slow=correction_on_slow_timescale)
thdforce.addForce(itype=3, filename="table_tf.xvg", type=typeCG)
integrator.addExtension(thdforce)

# distribute atoms and CG molecules according to AdResS domain decomposition, place CG molecules in the center of mass
espressopp.tools.AdressDecomp(system, integrator)

###########################
#  7. print system information  #
###########################

# print simulation parameters
print ''
print 'this is Force-based AdResS!'
print 'number of atomistic particles =', num_particles_AT
print 'number of coarse-grained particles =', num_particles_CG
print 'size =', size
print 'density =', density
print 'size of high resolution atomistic region =', ex_size
print 'size of hybrid region =', hy_size
print 'atomistic interaction cutoff =', interaction_cutoff_at
print 'coarse-grained interaction cutoff =', interaction_cutoff_cg
print 'atomistic interaction cutoff in verletlist =', interaction_cutoff_verletlist_at
print 'skin =', skin
print 'short timestep =', timestep
print 'long timestep =', longstep
print 'total run time =', totaltime
print 'total steps =', steps
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print ''

##################
#  8. run simulation  #
##################

# temperature analysis, timer, logfile
# note that the temperature will differ from target temperature since degrees of freedom are removed in the coarse-grained area, which is not taken into account here
temperature = espressopp.analysis.Temperature(system)
start_time = time.clock()
outfile = open("esp.dat", "w")

# output format
print 'Starting the integration loop...'
print ''
print "time, temperature, E_bonds, E_angles, E_nonbonded, E_kinetic, E_total"
fmt = '%5.5f %15.8g %15.10g %15.10g %15.10g %15.10g %15.10g\n'

# initial state
temp = temperature.compute()
E_bonds = 0
E_angles = 0
for bd in bondedinteractions.values(): E_bonds += bd.computeEnergy()
for ang in angleinteractions.values(): E_angles += ang.computeEnergy()
E_nonbonded = interNBat.computeEnergy() + interNBcg.computeEnergy()
E_kinetic = 0.5 * temp * (3 * num_particles_AT)
E_total = E_kinetic + E_bonds + E_angles + E_nonbonded
outfile.write(fmt%(0, temp, E_bonds, E_angles, E_nonbonded, E_kinetic, E_total))
print (fmt%(0, temp, E_bonds, E_angles, E_nonbonded, E_kinetic, E_total))

# run integration
for i in range(check):
    integrator.run(steps/check)
    temp = temperature.compute()
    E_bonds = 0
    E_angles = 0
    for bd in bondedinteractions.values(): E_bonds += bd.computeEnergy()
    for ang in angleinteractions.values(): E_angles += ang.computeEnergy()
    E_nonbonded = interNBat.computeEnergy() + interNBcg.computeEnergy()
    E_kinetic = 0.5 * temp * (3 * num_particles_AT)
    E_total = E_kinetic + E_bonds + E_angles + E_nonbonded
    outfile.write(fmt%((i+1)*steps*multistep/check*timestep, temp, E_bonds, E_angles, E_nonbonded, E_kinetic, E_total))
    print (fmt%((i+1)*steps*multistep/check*timestep, temp, E_bonds, E_angles, E_nonbonded, E_kinetic, E_total))

###########
#  9. Done  #
###########

# close log file
outfile.close()

# write the output configuration
espressopp.tools.writexyz("output.xyz", system)

# time Information
end_time = time.clock()
print 'Successfully finished simulation.'
print 'Run time = %.1f seconds' % (end_time - start_time)
