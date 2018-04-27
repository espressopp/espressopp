#!/usr/bin/env python
#  Copyright (C) 2017,2018
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

################################################################################
#                                                                                                                                            #
#  ESPResSo++ Python script for a quantum/classical path integral-based AdResS simulation  #
#                                                                                                                                            #
################################################################################

import time as pytime
import espressopp
import mpi4py.MPI as MPI
from espressopp import Real3D
from espressopp.tools import decomp

# Performs a path integral-based quantum/classical adaptive resolution simulation
# following the methodology proposed in J. Chem. Phys 147, 244104 (2017)
# The example system is liquid water, modeled using the force field proposed
# in J. Chem. Theory Comput. 10, 816 (2014). In the classical region a WCA potential
# between the oxygens is used.The simulation setup is similar to those used in the JCP paper.
print 'Performing an ESPResSo++ path integral-based quantum/classical adaptive resolution simulation.\n'

########################################################
#  1. specification of the system setup and simulation parameters  #
########################################################

print 'Specifying simulation parameters...'
steps =1000 # integration steps (outer steps of the multiple time stepping integrator)
intervals = 100 # intervals
timestep_short = 0.002/40.0 # the shortest timestep (interaction between Trotter beads on the ring polymers)
multiplier_short_to_medium = 10 # timestep_short * multiplier_short_to_medium gives the medium timestep (bonded interactions between atoms)
multiplier_medium_to_long = 4 # timestep_short * multiplier_short_to_medium * multiplier_medium_to_long gives the long timestep (non-bonded interactions between atoms)

interaction_cutoff = 0.84 # interaction cutoff for verletlist construction. Note: the verletlist is constucted using the ring polymers' centroids while the interaction takes place between the beads. Hence, this needs to be slightly larger (at least twice the ring polymers' radius of gyration) than the true potential cutoff
potential_cutoff = 0.78 # potential cutoff
skin = 0.1 # Verlet list skin

gamma = 2.0 # thermostat gamma
temp = 2.50266751 # thermostat temperature in reduced units

ex_size = 1.0 # redius of path integral region
hy_size = 1.5 # width of hybrid region

nTrotter = 32 # Trotter number
clmassmultiplier = 100.0 # classical mass multiplier with respect to real mass
constkinmass = False # speedup in classical region for interatomic interactions (centroids instead of all beads)?
PILE = True # Path Integral Langevin Equation thermostating (J. Chem. Phys. 133, 124104 (2010))? Only applicable when using real kinetic masses
PILElambda = 0.5 # PILE lambda parameter
realkinmass = True # real kinetic masses (for TRPMD)?
centroidthermostat = True # thermostat also the centroid (for CMD)?
CMDparameter = 1.0 # adiabatic decoupling gamma^2 for CMD
KTI = False # Kirkwood Thermodynamic Integration flag

speedupInterAtom = True # speedup in classical region for interatomic interactions (use centroids instead of all beads)?
speedupFreezeRings = False # speedup in classical region by freezing rings (skipping internal ring motion integration)? Note that, if set to True, the temperature will not be correct anymore, because freezing the rings in the CL region effectively reduces the number of degrees of freedom in the system

######################
#  2. read in coordinates  #
######################

print 'Reading in coordinates...'
# read equilibrated configuration file
pid, types, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz("input.xyz")

# Make masses
masses = []
for item in types:
    if item == 1:
        masses.append(15.9994)
    else:
        masses.append(1.008)

# Tables for Free Energy Correction and Thermodynamic Force on hydrogen and oxygen atoms
tabFEC_H = "FEC_H.dat"
tabFEC_O = "FEC_O.dat"
tabTHDF_H = "ThdForce_H.dat"
tabTHDF_O = "ThdForce_O.dat"

# Tables for potentials (see J. Chem. Theory Comput. 10, 816 (2014))
tabAngle = "POTS/tableESP_angle.dat"
tabBondHH = "POTS/tableESP_bondHH.dat"
tabBondOH = "POTS/tableESP_bondOH.dat"
tabHW_HW = "POTS/tableESP_HW_HW.dat"
tabHW_OW = "POTS/tableESP_HW_OW.dat"
tabOW_OW = "POTS/tableESP_OW_OW.dat"

num_Trotter_beads = len(x) # total number of Trotter beads in the system
num_atoms = len(x)/nTrotter # total number of atoms in the system
size = (Lx, Ly, Lz) # size

#####################
#  3. set up the system  #
#####################

print 'Setting up system...'
# System, boundary conditions, skin, communicator, node & cell grids
system = espressopp.System()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin
comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, interaction_cutoff, skin)

# Random number generator
rng = espressopp.esutil.RNG()
rng.seed(42)
system.rng = rng

# AdResS domain decomposition
system.storage = espressopp.storage.DomainDecompositionAdress(system, nodeGrid, cellGrid)

##########################
#  4. add particles to system  #
##########################

print 'Adding particles and tuples...'
props = ['id', 'pos', 'v', 'f', 'pib', 'type', 'mass', 'adrat']
allParticlesAT = []
allParticles = []
tuples = []

# prepare trotter beads (add here these particles just temporarily)
for pid_trotter in range(num_Trotter_beads):
    allParticlesAT.append([pid_trotter + 1,
                        Real3D(x[pid_trotter], y[pid_trotter], z[pid_trotter]), # position
                        Real3D(vx[pid_trotter], vy[pid_trotter], vz[pid_trotter]), # velocity
                        Real3D(0, 0, 0), # force
                        pid_trotter%nTrotter + 1, types[pid_trotter], masses[pid_trotter], 1]) # pib, type, mass, is AT particle

# create atoms
for pid_atom in range(num_atoms):

    # Preparation of tuples (tuples define, which atoms/trotter beads belong to which CG molecules/atoms)
    tmptuple = [pid_atom+num_Trotter_beads+1]
    for pid_trotter in range(nTrotter):
        pid = pid_atom*nTrotter+pid_trotter
        tmptuple.append((allParticlesAT[pid])[0])
    firstParticleId=tmptuple[1]
    cmp=allParticlesAT[firstParticleId-1][1]
    cmv=allParticlesAT[firstParticleId-1][2]

    # append atomistic particles
    allParticles.append([pid_atom+num_Trotter_beads+1,
                         Real3D(cmp[0], cmp[1], cmp[2]), # pos
                         Real3D(cmv[0], cmv[1], cmv[2]), # vel
                         Real3D(0, 0, 0), # force
                         0, types[pid_atom*nTrotter], masses[pid_atom*nTrotter], 0]) # pib, type, mass, is not AT particle
    # append Trotter beads
    for pid_trotter in range(nTrotter):
        pid = pid_atom*nTrotter+pid_trotter
        allParticles.append([(allParticlesAT[pid])[0],
                            (allParticlesAT[pid])[1], # pos
                            (allParticlesAT[pid])[2], # vel
                            (allParticlesAT[pid])[3], # force
                            (allParticlesAT[pid])[4], # pib
                            (allParticlesAT[pid])[5], # type
                            (allParticlesAT[pid])[6], # mass
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
#  4. set up structure, interactions, and force field  #
###########################################

print 'Setting up interactions and force field...'
# create bond lists between atoms
bondsOH = []
bondsHH = []
for part in range(num_atoms/3):
    bondsOH.append((num_Trotter_beads + 1 + 3*part, num_Trotter_beads + 1 + 3*part+1))
    bondsOH.append((num_Trotter_beads + 1 + 3*part, num_Trotter_beads + 1 + 3*part+2))
    bondsHH.append((num_Trotter_beads + 1 + 3*part+1, num_Trotter_beads + 1 + 3*part+2))

# add bonds between atoms
fplOH = espressopp.FixedPairList(system.storage)
fplHH = espressopp.FixedPairList(system.storage)
fplOH.addBonds(bondsOH)
fplHH.addBonds(bondsHH)

# create the adaptive resolution verlet list
vl = espressopp.VerletListAdress(system, cutoff=interaction_cutoff, adrcut=interaction_cutoff,
                                dEx=ex_size, dHy=hy_size,
                                adrCenter=[Lx/2, Ly/2, Lz/2], exclusionlist=bondsOH+bondsHH)

# create angle list between atoms
angles = []
for part in range(num_atoms/3):
    angles.append((num_Trotter_beads + 1 + 3*part+1, num_Trotter_beads + 1 + 3*part, num_Trotter_beads + 1 + 3*part+2))

# add angles between atoms
ftl = espressopp.FixedTripleList(system.storage)
ftl.addTriples(angles)

# non-bonded potentials
# inteaction first (we need specific PI AdResS interaction type)
interNB = espressopp.interaction.VerletListPIadressTabulatedLJ(vl, ftpl, nTrotter, speedupInterAtom)
# QM potential
potOOqm = espressopp.interaction.Tabulated(itype=3, filename=tabOW_OW, cutoff=potential_cutoff)
potHOqm = espressopp.interaction.Tabulated(itype=3, filename=tabHW_OW, cutoff=potential_cutoff)
potHHqm = espressopp.interaction.Tabulated(itype=3, filename=tabHW_HW, cutoff=potential_cutoff)
interNB.setPotentialQM(type1=1, type2=1, potential=potOOqm)
interNB.setPotentialQM(type1=1, type2=0, potential=potHOqm)
interNB.setPotentialQM(type1=0, type2=0, potential=potHHqm)
# WCA potential
potOOcl = espressopp.interaction.LennardJones(epsilon=temp, sigma=0.25, shift='auto', cutoff=1.122462048309373*0.25)
interNB.setPotentialCL(type1=1, type2=1, potential=potOOcl)
# add interaction to system
system.addInteraction(interNB)

# bonded potentials
potBondHH = espressopp.interaction.Tabulated(itype=3, filename=tabBondHH)
potBondOH = espressopp.interaction.Tabulated(itype=3, filename=tabBondOH)
interBondedHH = espressopp.interaction.FixedPairListPIadressTabulated(system, fplHH, ftpl, potBondHH, nTrotter, speedupInterAtom)
interBondedOH = espressopp.interaction.FixedPairListPIadressTabulated(system, fplOH, ftpl, potBondOH, nTrotter, speedupInterAtom)
system.addInteraction(interBondedHH)
system.addInteraction(interBondedOH)

# angle potentials
potAngle = espressopp.interaction.TabulatedAngular(itype=3, filename=tabAngle)
interAngle = espressopp.interaction.FixedTripleListPIadressTabulatedAngular(system, ftl, ftpl, potAngle, nTrotter, speedupInterAtom)
system.addInteraction(interAngle)

#########################################
#  5. set up integration scheme and corrections  #
#########################################

print 'Setting up integration scheme (path integral-based adaptive resolution multiple time stepping integrator)...'
# path integral-based adaptive resolution multiple time stepping integrator
integrator = espressopp.integrator.PIAdressIntegrator(system=system, verletlist=vl, timestep=timestep_short, sSteps=multiplier_short_to_medium, mSteps=multiplier_medium_to_long, nTrotter=nTrotter, realKinMass=realkinmass, constKinMass=constkinmass, temperature=temp, gamma=gamma, centroidThermostat=centroidthermostat, CMDparameter=CMDparameter, PILE=PILE, PILElambda=PILElambda, CLmassmultiplier=clmassmultiplier, speedup=speedupFreezeRings, KTI=KTI)

# add Free Energy Correction for oxygen and hydrogen atoms
fec = espressopp.integrator.FreeEnergyCompensation(system, center=[Lx/2, Ly/2, Lz/2], ntrotter=nTrotter)
fec.addForce(itype=3, filename=tabFEC_O, type=1)
fec.addForce(itype=3, filename=tabFEC_H, type=0)
integrator.addExtension(fec)

# add Thermodynamic Force for oxygen and hydrogen atoms
thdf = espressopp.integrator.TDforce(system, vl)
thdf.addForce(itype=3, filename=tabTHDF_O, type=1)
thdf.addForce(itype=3, filename=tabTHDF_H, type=0)
integrator.addExtension(thdf)

# distribute path integral beads and atoms according to AdResS domain decomposition, place AT/CG particles into the ring polymers' centers of mass
espressopp.tools.AdressDecomp(system, integrator)

##########################
#  6. set up analysis routines  #
##########################

print 'Setting up analysis routines...'
# radius of gyration profiles
gyration_array_total_H = []
gyration_array_total_O = []
gyrationprofile = espressopp.analysis.RadGyrXProfilePI(system)
gyrationprofilegrid = 40
gyrationAdds_H = [0 for i in range(gyrationprofilegrid)]
gyrationAdds_O = [0 for i in range(gyrationprofilegrid)]

# OO rdf
rdf_array_total_OO = []
Adds_OO = 0.0
rdf_OO = espressopp.analysis.RDFatomistic(system, 1, 1, 0.75)

# OH rdf
rdf_array_total_OH = []
Adds_OH = 0.0
rdf_OH = espressopp.analysis.RDFatomistic(system, 1, 0, 0.75)

# HH rdf
rdf_array_total_HH = []
Adds_HH = 0.0
rdf_HH = espressopp.analysis.RDFatomistic(system, 0, 0, 0.75)

# grids and formatting for analysis routines
rdfgrid = 400
fmt_rdf = ' %12.8f %12.8f\n'
dr_rdf = Ly / (2.0*float(rdfgrid))
fmt_gyr = ' %12.8f %12.8f\n'
dr_gyr = Lx / float(gyrationprofilegrid)

# output log
outfile = open("esp.dat", "w")

###########################
#  7. print system information  #
###########################

print ''
print "System setup done, information:"
print ''
print 'PI-AdResS Center =', [Lx/2, Ly/2, Lz/2]
print 'Size of high resolution full path integral region', ex_size
print 'Size of hybrid region', hy_size
print 'Trotter number =', integrator.getNtrotter()
print 'Total number of Trotter beads =', num_Trotter_beads
print 'Total number of atoms =', num_atoms
print 'Atomistic density = %.4f' % (num_atoms / (Lx * Ly * Lz))
print 'Interaction cutoff =', interaction_cutoff
print 'Potential cutoff =', potential_cutoff
print 'Skin =', system.skin
print 'Short timestep =', integrator.getTimeStep()
print 'Medium timestep =', integrator.getTimeStep() * integrator.getsStep()
print 'Long timestep =', integrator.getTimeStep() * integrator.getmStep() * integrator.getsStep()
print 'Outer steps =', steps
print 'Intervals =', intervals
print 'NodeGrid = %s' % (nodeGrid,)
print 'CellGrid = %s' % (cellGrid,)
print 'Temperature =', integrator.getTemperature()
print 'Gamma =', integrator.getGamma()
print 'Classical Mass Multiplier =', integrator.getClmassmultiplier()
print 'Constant kinetic mass?', integrator.getConstKinMass()
print 'Using (adaptive or constant) real kinetic masses?', integrator.getRealKinMass()
print 'Path Integration Langevin Equation thermostating?', integrator.getPILE()
print 'Path Integration Langevin Equation thermostating lambda =', integrator.getPILElambda()
print 'Thermostating the centroid?', integrator.getCentroidThermostat()
print 'CMD adiabadicity parameter =', integrator.getCMDparameter()
print 'Running Kirkwood Thermodynamic Integration?', integrator.getKTI()
print 'Using centers of mass in classical region for force calculations?', speedupInterAtom
print 'Freezing internal ring vibrations in classical region?', integrator.getSpeedup()
print ''

##################
#  8. run simulation  #
##################

# timer, steps
nsteps = steps / intervals
start_time = pytime.clock()

# output format for screen and file
print 'Starting the integration loop...'
print ''
print 'step, time (ps), temperature, E_bonds, E_angles, E_ringpolymer, E_nonbonded, E_kin, E_correction, E_total'
fmt = '%8d %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g\n'

# initial configuration analysis
Eb = interBondedOH.computeEnergy() + interBondedHH.computeEnergy()
EAng = interAngle.computeEnergy()
ELj= interNB.computeEnergy()
Ek = integrator.computeKineticEnergy()
T = Ek *2.0 / (0.00831451 * 3.0 * num_Trotter_beads)
EPI = integrator.computeRingEnergy()
Ecorr = fec.computeCompEnergy() + thdf.computeTDEnergy()
Etotal = Ek+Eb+EAng+ELj+EPI+Ecorr
outfile.write(fmt%(0, 0, T, Eb, EAng, EPI, ELj, Ek, Ecorr, Etotal))
print ((fmt%(0, 0, T, Eb, EAng, EPI, ELj, Ek, Ecorr, Etotal)))

# integration and on the fly analysis
for s in range(1, intervals + 1):
    integrator.run(nsteps)
    step = nsteps * s
    time = step * timestep_short * multiplier_medium_to_long * multiplier_short_to_medium
    Eb = interBondedOH.computeEnergy() + interBondedHH.computeEnergy()
    EAng = interAngle.computeEnergy()
    ELj= interNB.computeEnergy()
    Ek = integrator.computeKineticEnergy()
    T = Ek *2.0 / (0.00831451 * 3.0 * num_Trotter_beads)
    EPI = integrator.computeRingEnergy()
    Ecorr = fec.computeCompEnergy() + thdf.computeTDEnergy()
    Etotal = Ek+Eb+EAng+ELj+EPI+Ecorr
    outfile.write(fmt%(step, time, T, Eb, EAng, EPI, ELj, Ek, Ecorr, Etotal))
    print (fmt%(step, time, T, Eb, EAng, EPI, ELj, Ek, Ecorr, Etotal))

    if s%10==0:
        rdf_array_OO = rdf_OO.computePathIntegral(rdfgrid)
        for i in range(len(rdf_array_OO)):
                if(i>=len(rdf_array_total_OO)):
                        rdf_array_total_OO.append(rdf_array_OO[i])
                else:
                        rdf_array_total_OO[i] += rdf_array_OO[i]
        Adds_OO += 1.0

        rdf_array_OH = rdf_OH.computePathIntegral(rdfgrid)
        for i in range(len(rdf_array_OH)):
                if(i>=len(rdf_array_total_OH)):
                        rdf_array_total_OH.append(rdf_array_OH[i])
                else:
                        rdf_array_total_OH[i] += rdf_array_OH[i]
        Adds_OH += 1.0

        rdf_array_HH = rdf_HH.computePathIntegral(rdfgrid)
        for i in range(len(rdf_array_HH)):
                if(i>=len(rdf_array_total_HH)):
                        rdf_array_total_HH.append(rdf_array_HH[i])
                else:
                        rdf_array_total_HH[i] += rdf_array_HH[i]
        Adds_HH += 1.0

    gyration_array_H = gyrationprofile.compute(gyrationprofilegrid, nTrotter, 0)
    for i in range(len(gyration_array_H)):
                if(i>=len(gyration_array_total_H)):
                        gyration_array_total_H.append(gyration_array_H[i])
                else:
                        gyration_array_total_H[i] += gyration_array_H[i]
                if(gyration_array_H[i] != 0.0):
                    gyrationAdds_H[i] += 1.0

    gyration_array_O = gyrationprofile.compute(gyrationprofilegrid, nTrotter, 1)
    for i in range(len(gyration_array_O)):
                if(i>=len(gyration_array_total_O)):
                        gyration_array_total_O.append(gyration_array_O[i])
                else:
                        gyration_array_total_O[i] += gyration_array_O[i]
                if(gyration_array_O[i] != 0.0):
                    gyrationAdds_O[i] += 1.0

# close output file
outfile.close()

###################################################
#  8. postprocess gathered analysis outputs and print to file  #
###################################################

# print O-O rdf to file
for i in range(len(rdf_array_total_OO)):
  rdf_array_total_OO[i] /= Adds_OO
rdf_OO_file = open('rdf_profile_OO.dat', 'w')
for i in range(len(rdf_array_total_OO)):
  rdf_OO_file.write(fmt_rdf % ( (i+0.5)*dr_rdf, rdf_array_total_OO[i] ))
rdf_OO_file.close()

# print O-H rdf to file
for i in range(len(rdf_array_total_OH)):
  rdf_array_total_OH[i] /= Adds_OH
rdf_OH_file = open('rdf_profile_OH.dat', 'w')
for i in range(len(rdf_array_total_OH)):
  rdf_OH_file.write(fmt_rdf % ( (i+0.5)*dr_rdf, rdf_array_total_OH[i] ))
rdf_OH_file.close()

# print H-H rdf to file
for i in range(len(rdf_array_total_HH)):
  rdf_array_total_HH[i] /= Adds_HH
rdf_HH_file = open('rdf_profile_HH.dat', 'w')
for i in range(len(rdf_array_total_HH)):
  rdf_HH_file.write(fmt_rdf % ( (i+0.5)*dr_rdf, rdf_array_total_HH[i] ))
rdf_HH_file.close()

# print hydrogen radius of gyration profile to file
for i in range(len(gyration_array_total_H)):
  if(gyrationAdds_H[i] > 0.0):
    gyration_array_total_H[i] /= gyrationAdds_H[i]
rgyr_H_file = open('radgyr_profile_H.dat', 'w')
for i in range(len(gyration_array_total_H)):
  rgyr_H_file.write(fmt_gyr % ( (i+0.5)*dr_gyr, gyration_array_total_H[i] ))
rgyr_H_file.close()

# print oxygen radius of gyration profile to file
for i in range(len(gyration_array_total_O)):
  if(gyrationAdds_O[i] > 0.0):
    gyration_array_total_O[i] /= gyrationAdds_O[i]
rgyr_O_file = open('radgyr_profile_O.dat', 'w')
for i in range(len(gyration_array_total_O)):
  rgyr_O_file.write(fmt_gyr % ( (i+0.5)*dr_gyr, gyration_array_total_O[i] ))
rgyr_O_file.close()

###########
#  9. Done  #
###########

end_time = pytime.clock()
print 'Successfully finished simulation.'
print 'Run time = %.1f seconds' % (end_time - start_time)
