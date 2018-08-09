#!/usr/bin/env python2
#  Copyright (C) 2012-2017(H)
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
#  ESPResSo++ Python script for tabulated GROMACS simulation              #
#                                                                         #
###########################################################################

import sys
import time
import espressopp
import mpi4py.MPI as MPI
import logging
from espressopp import Real3D, Int3D
from espressopp.tools import gromacs
from espressopp.tools import decomp
from espressopp.tools import timers


# simulation parameters (nvt = False is nve)
steps = 1 
check = steps
rc    = 1.1  # Verlet list cutoff
skin  = 0.3
timestep = 0.0004

# GROMACS tabulated potentials files
tabAAg = "table_A_A.xvg"    # non-bonded
tabABg = "table_A_B.xvg"
tabBBg = "table_B_B.xvg"

taba1g = "table_a1.xvg"     # angles
taba2g = "table_a2.xvg"
taba3g = "table_a3.xvg"  
taba4g = "table_a4.xvg"   
taba5g = "table_a5.xvg"

tabb0g  = "table_b0.xvg"    # bonds
tabb1g  = "table_b1.xvg"
tabb5g  = "table_b5.xvg"
tabb12g = "table_b12.xvg"
tabb17g = "table_b17.xvg"
#tabb18g = "table_b18.xvg"
tabb19g = "table_b19.xvg"
tabb40g = "table_b40.xvg"

tabd6g = "table_d6.xvg"     # dihedrals
tabd7g = "table_d7.xvg"
tabd8g = "table_d8.xvg"
tabd9g = "table_d9.xvg"

spline = 2                  # spline interpolation type (1, 2, 3)

tabfilesnb = [tabAAg, tabABg, tabBBg]
#tabfiles2b = [tabb0g, tabb1g, tabb5g, tabb12g, tabb17g, tabb19g, tabb40g]
#tabfiles3b = [taba1g, taba2g, taba3g, taba4g, taba5g]
#tabfiles4b = [tabd6g, tabd7g, tabd8g, tabd9g]


# parameters to convert GROMACS tabulated potential file
sigma = 1.0
epsilon = 1.0
c6 = 1.0
c12 = 1.0

# GROMACS setup files
grofile = "confout.gro"
topfile = "topol.top"

# Generate dictonary of potentials for all particle types (relies on table file names).
# Input is list of gromacs tabulated non-bonded potentials file names
# Output is something like {"A_A":potAA, "A_B":potAB, "B_B":potBB}
def genTabPotentials(tabfilesnb):
    potentials = {}
    for fg in tabfilesnb:
        fe = fg.split(".")[0]+".tab" # name of espressopp file
        gromacs.convertTable(fg, fe, sigma, epsilon, c6, c12)
        pot = espressopp.interaction.Tabulated(itype=spline, filename=fe, cutoff=rc)
        t1, t2 = fg[6], fg[8] # type 1, type 2
        potentials.update({t1+"_"+t2: pot})
    return potentials

# define types of particles (used for non-bonded interactions)
particleTypes = {"A":["A1m", "A2m", "A1r", "A2r"],"B":["B1u","B2u", "B1d", "B2d"]}
potentials = genTabPotentials(tabfilesnb)


######################################################################
##  IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE  ##
######################################################################
defaults, types, atomtypes, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, dihedraltypes, dihedraltypeparams, exclusions, x, y, z, vx, vy, vz, resname, resid, Lx, Ly, Lz =gromacs.read(grofile, topfile)
#defaults, types, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, dihedraltypes, dihedraltypeparams, exclusions, x, y, z, vx, vy, vz, Lx, Ly, Lz = gromacs.read(grofile, topfile)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

sys.stdout.write('Setting up simulation ...\n')
system = espressopp.System()
system.rng = espressopp.esutil.RNG()
system.bc = espressopp.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size,size,rc,skin)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# add particles to the system and then decompose
props = ['id', 'pos', 'v', 'type']
allParticles = []
for pid in range(num_particles):
    part = [pid + 1, Real3D(x[pid], y[pid], z[pid]),
            Real3D(vx[pid], vy[pid], vz[pid]), types[pid]]
    allParticles.append(part)
system.storage.addParticles(allParticles, *props)    
system.storage.decompose()



# Tabulated Verlet list for non-bonded interactions
vl = espressopp.VerletList(system, cutoff = rc + system.skin)
internb = espressopp.interaction.VerletListTabulated(vl)
gromacs.setTabulatedInteractions(potentials, particleTypes, system, internb)

vl.exclude(exclusions)


# bonded 2-body interactions
bondedinteractions=gromacs.setBondedInteractions(system, bondtypes, bondtypeparams)
# alternative, manual way for setting up bonded interactions
#for k, v in bonds.iteritems(): # k is number of potential table, v is bondlist
#    fpl = espressopp.FixedPairList(system.storage)
#    fpl.addBonds(v)
#    fg = "table_b"+k+".xvg"
#    fe = fg.split(".")[0]+".tab" # name of espressopp file
#    gromacs.convertTable(fg, fe, sigma, epsilon, c6, c12)
#    potTab = espressopp.interaction.Tabulated(itype=spline, filename=fe)
#    interb = espressopp.interaction.FixedPairListTabulated(system, fpl, potTab)
#    system.addInteraction(interb)


# bonded 3-body interactions
angleinteractions=gromacs.setAngleInteractions(system, angletypes, angletypeparams)
# alternative, manual way for setting up angular interactions
#for k, v in angles.iteritems(): # k is number of potential table, v is anglelist
#    ftl = espressopp.FixedTripleList(system.storage)
#    ftl.addTriples(v)
#    fg = "table_a"+k+".xvg"
#    fe = fg.split(".")[0]+".tab" # name of espressopp file
#    gromacs.convertTable(fg, fe, sigma, epsilon, c6, c12)
#    potTab = espressopp.interaction.TabulatedAngular(itype=spline, filename=fe)
#    intera = espressopp.interaction.FixedTripleListTabulatedAngular(system, ftl, potTab)
#    system.addInteraction(intera)
    

# bonded 4-body interactions
dihedralinteractions=gromacs.setDihedralInteractions(system, dihedraltypes, dihedraltypeparams)
# alternative, manual way for setting up dihedral interactions:
#for k, v in dihedrals.iteritems(): # k is number of potential table, v is anglelist
#    fql = espressopp.FixedQuadrupleList(system.storage)
#    fql.addQuadruples(v)
#    fg = "table_d"+k+".xvg"
#    fe = fg.split(".")[0]+".tab" # name of espressopp file
#    gromacs.convertTable(fg, fe, sigma, epsilon, c6, c12)
#    potTab = espressopp.interaction.TabulatedDihedral(itype=spline, filename=fe)
#    interd = espressopp.interaction.FixedQuadrupleListTabulatedDihedral(system, fql, potTab)
#    system.addInteraction(interd)



#exit()



# langevin thermostat
langevin = espressopp.integrator.LangevinThermostat(system)
langevin.gamma = 1.0
langevin.temperature = 1.0
integrator = espressopp.integrator.VelocityVerlet(system)
integrator.addExtension(langevin)
integrator.dt = timestep


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


#exit()



# analysis
configurations = espressopp.analysis.Configurations(system)
configurations.gather()
temperature = espressopp.analysis.Temperature(system)
pressure = espressopp.analysis.Pressure(system)
pressureTensor = espressopp.analysis.PressureTensor(system)

fmt = '%5d %8.4f %11.4f %11.4f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f\n'

T = temperature.compute()
#P = pressure.compute()
P = 0
#Pij = pressureTensor.compute()
Pij = [0,0,0,0,0,0]
Ek = 0.5 * T * (3 * num_particles)
Ep = internb.computeEnergy()
Eb, Ea, Ed=0,0,0
for bd in bondedinteractions.values():Eb+=bd.computeEnergy()
for ang in angleinteractions.values(): Ea+=ang.computeEnergy()
for dih in dihedralinteractions.values(): Ed+=dih.computeEnergy()

Etotal = Ek + Ep + Eb + Ea + Ed
sys.stdout.write(' step     T          P          Pxy        etotal      ekinetic      epair         ebond       eangle       edihedral\n')
sys.stdout.write(fmt % (0, T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea, Ed))

start_time = time.clock()

espressopp.tools.pdb.pdbwrite("traj.pdb", system, append=False)

for i in range(check):

    integrator.run(steps/check) # print out every steps/check steps
    T = temperature.compute()
    P = pressure.compute()
    Pij = pressureTensor.compute()
    Ek = 0.5 * T * (3 * num_particles)
    Ep = internb.computeEnergy()
    
    Eb, Ea, Ed=0,0,0
    for bd in bondedinteractions.values():Eb+=bd.computeEnergy()
    for ang in angleinteractions.values(): Ea+=ang.computeEnergy()
    for dih in dihedralinteractions.values(): Ed+=dih.computeEnergy()
    Etotal = Ek + Ep + Eb + Ea + Ed
    sys.stdout.write(fmt % ((i+1)*(steps/check), T, P, Pij[3], Etotal, Ek, Ep, Eb, Ea, Ed))
    #espressopp.tools.pdb.pdbwrite("traj.pdb", system, append=True)
    #sys.stdout.write('\n')

# print timings and neighbor list information
end_time = time.clock()
timers.show(integrator.getTimers(), precision=2)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))




