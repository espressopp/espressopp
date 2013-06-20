#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################################################
#                                                                         #
#  ESPResSo++ Python script for tabulated GROMACS simulation              #
#                                                                         #
###########################################################################

import sys
import time
import espresso
import MPI
import logging
from espresso import Real3D, Int3D
from espresso.tools.convert import gromacs
from espresso.tools import decomp
from espresso.tools import timers


# simulation parameters (nvt = False is nve)
steps = 1000
check = 6   # how many times to display energies during run
rc   = 2.5  # Verlet list cutoff
rcaa = 2.5  # cutoff A-A
rcab = 2.0  # cutoff A-B
rcbb = 2.0  # cutoff B-B
skin = 0.3
timestep = 0.001

# GROMACS tabulated potentials files
tabAAg = "table_A_A.xvg"        # non-bonded
tabABg = "table_A_B.xvg"
tabBBg = "table_B_B.xvg"

tab2bg = "table_b1.xvg"         # 2-body bonded
tab3bg = "table_a1.xvg"         # 3-body bonded
spline = 2                      # spline interpolation type (1, 2, 3)



# parameters to convert GROMACS tabulated potential file
tabAA = "table_A_A.tab"         # output file names
tabAB = "table_A_B.tab"
tabBB = "table_B_B.tab"
tab2b = "table_b1.tab" 
tab3b = "table_a1.tab"
sigma = 1.0
epsilon = 1.0
c6 = 1.0
c12 = 1.0

# GROMACS setup files
grofile = "conf.gro"
topfile = "topol.top"




######################################################################
##  IT SHOULD BE UNNECESSARY TO MAKE MODIFICATIONS BELOW THIS LINE  ##
######################################################################
defaults, types, masses, charges, atomtypeparameters, bondtypes, bondtypeparams, angletypes, angletypeparams, exclusions, x, y, z, Lx, Ly, Lz= gromacs.read(grofile, topfile)
num_particles = len(x)
density = num_particles / (Lx * Ly * Lz)
size = (Lx, Ly, Lz)

sys.stdout.write('Setting up simulation ...\n')
system = espresso.System()
system.rng = espresso.esutil.RNG()
system.bc = espresso.bc.OrthorhombicBC(system.rng, size)
system.skin = skin

comm = MPI.COMM_WORLD
nodeGrid = decomp.nodeGrid(comm.size)
cellGrid = decomp.cellGrid(size, nodeGrid, rc, skin)
system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)


# add particles to the system and then decompose
for pid in range(num_particles):
    #system.storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
    system.storage.addParticles([[pid + 1, Real3D(x[pid], y[pid], z[pid]), types[pid]]], "id", "pos", "type")
system.storage.decompose()


# convert gromacs tabulated files to espresso++ format
gromacs.convertTable(tabAAg, tabAA, sigma, epsilon, c6, c12)
gromacs.convertTable(tabABg, tabAB, sigma, epsilon, c6, c12)
gromacs.convertTable(tabBBg, tabBB, sigma, epsilon, c6, c12)
gromacs.convertTable(tab2bg, tab2b, sigma, epsilon, c6, c12)
gromacs.convertTable(tab3bg, tab3b, sigma, epsilon, c6, c12)



# non-bonded interactions, B is type 0, A is type 1
# Verlet list
vl = espresso.VerletList(system, cutoff = rc + system.skin)
# note: in the previous version of this example, exclusions were treated
# incorrectly. Here the nrexcl=3 parameter is taken into account 
# which excludes all neighbors up to 3 bonds away
vl.exclude(exclusions)

internb = espresso.interaction.VerletListTabulated(vl)
# A-A with Verlet list      
potTab = espresso.interaction.Tabulated(itype=spline, filename=tabAA, cutoff=rcaa)
internb.setPotential(type1 = 1, type2 = 1, potential = potTab)
# A-B with Verlet list      
potTab = espresso.interaction.Tabulated(itype=spline, filename=tabAB, cutoff=rcab)
internb.setPotential(type1 = 1, type2 = 0, potential = potTab)
# B-B with Verlet list      
potTab = espresso.interaction.Tabulated(itype=spline, filename=tabBB, cutoff=rcbb)
internb.setPotential(type1 = 0, type2 = 0, potential = potTab)
system.addInteraction(internb)



# 2-body bonded interactions
bondedinteractions=gromacs.setBondedInteractions(system, bondtypes, bondtypeparams)
#This could also be done manually:
#fpl = espresso.FixedPairList(system.storage)
#fpl.addBonds(bonds['1'])
#potTab = espresso.interaction.Tabulated(itype=spline, filename=tab2b)
#interb = espresso.interaction.FixedPairListTabulated(system, fpl, potTab)
#system.addInteraction(interb)

# 3-body bonded interactions
angleinteractions=gromacs.setAngleInteractions(system, angletypes, angletypeparams)

#This could also be done manually:
#ftl = espresso.FixedTripleList(system.storage)
#ftl.addTriples(angles['1'])
#potTab = espresso.interaction.TabulatedAngular(itype=spline, filename = tab3b)
#intera = espresso.interaction.FixedTripleListTabulatedAngular(system, ftl, potTab)
#system.addInteraction(intera)



# langevin thermostat
langevin = espresso.integrator.LangevinThermostat(system)
langevin.gamma = 1.0
langevin.temperature = 1.0
integrator = espresso.integrator.VelocityVerlet(system)
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
configurations = espresso.analysis.Configurations(system)
configurations.gather()
temperature = espresso.analysis.Temperature(system)
pressure = espresso.analysis.Pressure(system)
pressureTensor = espresso.analysis.PressureTensor(system)

fmt = '%5d %8.4f %11.4f %11.4f %12.3f %12.3f %12.3f %12.3f %12.3f\n'

T = temperature.compute()
P = pressure.compute()
Pij = pressureTensor.compute()
Ek = 0.5 * T * (3 * num_particles)
#Epaa = interaa.computeEnergy()
#Epab = interab.computeEnergy()
#Epbb = interbb.computeEnergy()
Epnb = internb.computeEnergy()
Eb = bondedinteractions[0].computeEnergy()
Ea = angleinteractions[0].computeEnergy()
Etotal = Ek + Epnb + Eb + Ea
sys.stdout.write(' step     T          P          Pxy        etotal      ekinetic      epair        ebond       eangle\n')
sys.stdout.write(fmt % (0, T, P, Pij[3], Etotal, Ek, Epnb, Eb, Ea))

start_time = time.clock()

for i in range(check):
    integrator.run(steps/check)
    T = temperature.compute()
    P = pressure.compute()
    Pij = pressureTensor.compute()
    Ek = 0.5 * T * (3 * num_particles)
    Epnb = internb.computeEnergy()
    Eb = bondedinteractions[0].computeEnergy()
    Ea = angleinteractions[0].computeEnergy()
    Etotal = Ek + Epnb + Eb + Ea
    sys.stdout.write(fmt % ((i+1)*(steps/check), T, P, Pij[3], Etotal, Ek, Epnb, Eb, Ea))
    #sys.stdout.write('\n')

# print timings and neighbor list information
end_time = time.clock()
timers.show(integrator.getTimers(), precision=2)
sys.stdout.write('Total # of neighbors = %d\n' % vl.totalSize())
sys.stdout.write('Ave neighs/atom = %.1f\n' % (vl.totalSize() / float(num_particles)))
sys.stdout.write('Neighbor list builds = %d\n' % vl.builds)
sys.stdout.write('Integration steps = %d\n' % integrator.step)
sys.stdout.write('CPU time = %.1f\n' % (end_time - start_time))




