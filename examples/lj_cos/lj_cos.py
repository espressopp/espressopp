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
# DEMONSTRATION OF THE LJ-COS POTENTIAL CONTROLLING SOLVENT QUALITY       #
# BY PARAMETER PHI (see M. Steinhauser PhD thesis for details)            #  
###########################################################################

import espressopp
from espressopp import Int3D
from espressopp import Real3D

# initial parameters of the simulation
num_chains     = 5
mon_per_chain  = 100
L              = 20
box            = (L, L, L)
bondlen        = 0.97
rc             = 1.5
skin           = 2.5
dt             = 0.00001
sigma          = 1.
temperature    = 1.
phi            = 0.

system         = espressopp.System()
system.rng     = espressopp.esutil.RNG()
system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
nodeGrid       = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size,box,rc,skin)
cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# set-up of the integrator and its timestep
integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = dt

# set-up of the thermostat
thermostat     = espressopp.integrator.LangevinThermostat(system)
thermostat.gamma  = 1.0
thermostat.temperature = temperature
integrator.addExtension(thermostat)

print 'timestep is ', integrator.dt
print 'gamma of the thermostat is ', thermostat.gamma
print 'temperature of the thermostat is ', thermostat.temperature

props    = ['id', 'type', 'mass', 'pos', 'v']
vel_zero = espressopp.Real3D(0.0, 0.0, 0.0)

bondlist = espressopp.FixedPairList(system.storage)
pid      = 1
ptype    = 0
mass     = 1.0
chain    = []

exclusionlist = [] # to separate bonded and non-bonded LJ potentials; used for equilibration
for i in range(num_chains):
    startpos = system.bc.getRandomPos()
    positions, bonds = espressopp.tools.topology.polymerRW(pid, startpos, mon_per_chain, bondlen)
    for k in range(mon_per_chain):
        part = [pid + k, ptype, mass, positions[k], vel_zero]
        chain.append(part)
        if (k < mon_per_chain - 1): 
            exclusionlist.append(bonds[k])
    pid += mon_per_chain
    system.storage.addParticles(chain, *props)
    system.storage.decompose()
    chain = []
    bondlist.addBonds(bonds)
    
system.storage.decompose()

# initial set-up of the Lennard-Jones cosine
vl = espressopp.VerletList(system, cutoff = rc)
vl.exclude(exclusionlist)

# non-bonded LJcos potential
nbLJcos = espressopp.interaction.LJcos(phi=phi)
nbLJcosInter = espressopp.interaction.VerletListLJcos(vl)
nbLJcosInter.setPotential(type1=0, type2=0, potential=nbLJcos)
system.addInteraction(nbLJcosInter)

# bonded LJcos potential
bondLJcos = espressopp.interaction.LJcos(phi=phi)
bondLJcos.sigma = sigma
bondLJcosInter = espressopp.interaction.FixedPairListLJcos(system, bondlist, bondLJcos)
system.addInteraction(bondLJcosInter)

# FENE-potential
potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, bondlist, potFENE)
system.addInteraction(interFENE)

def main ():
    warmUpFirst()
    warmUpSecond()

    print "Starting simulation" 
    for k in range(100):
        integrator.run(1000)
    
        findTemperature()
        findRee()

        espressopp.tools.analyse.info(system, integrator)

def warmUpFirst():
    # FIRST PHASE OF THE WARM UP
    # force-capping
    max_force = 10000.
    force_capping = espressopp.integrator.CapForce(system, max_force)
    integrator.addExtension(force_capping)

    print "First phase of the warm up. Sigma will be increased from 0. to 1.0 and timestep to 0.001"
    new_sigma = 0.
    for l in range(4):
        print "start increasing sigma from ", new_sigma
        for k in range(5):
            new_sigma += 0.05
            nbLJcos.sigma = new_sigma
            nbLJcosInter.setPotential(type1=0, type2=0, potential=nbLJcos)
            integrator.run(1000)
            espressopp.tools.analyse.info(system, integrator)
        print "increased sigma to ", new_sigma
        print "increasing timestep from ", integrator.dt
        dt_step = .20 * integrator.dt
        for k in range(10):
            integrator.dt += dt_step
            integrator.run(1000)
            espressopp.tools.analyse.info(system, integrator)
        print "increased timestep to ", integrator.dt

    force_capping.disconnect()
    print "switching off force capping"

    print "new_sigma is ", new_sigma
    # FINISHED THE FIRST PHASE OF THE WARM UP

def warmUpSecond():
    # SECOND PHASE OF THE WARM UP
    integrator.step = 0
    integrator.dt = 0.005
    print "Second phase of the warm up with a production timestep of", integrator.dt, ". Force capping is turned off."
    for k in range(100):
        integrator.run(10000)
        espressopp.tools.analyse.info(system, integrator)

    num_particles = num_chains * mon_per_chain

    integrator.step = 0
    # FINISHED THE SECOND PHASE OF THE WARM UP

def findTemperature():
    T   = espressopp.analysis.Temperature(system)

    # write temperature
    f_temp = open('temperature.dat', 'a')

    s = str(integrator.step)
    f_temp.write(s+'\t')

    currT = T.compute()
    s = str(currT)
    f_temp.write(s+'\n')

    f_temp.close()

def findRee():
    # write end-to-end distance and gyration radius
    f_stat = open('static_prop.dat', 'a')

    s = str(integrator.step)
    f_stat.write(s+'\t')

    # calculate end-to-end
    R2ee = 0.
    R2gyr = 0.

    global num_chains, mon_per_chain
    for i in range(num_chains):
        p1 = system.storage.getParticle(i*mon_per_chain+1)
        p2 = system.storage.getParticle((i+1)*mon_per_chain)
        c1 = system.bc.getUnfoldedPosition(p1.pos, p1.imageBox)
        c2 = system.bc.getUnfoldedPosition(p2.pos, p2.imageBox)
        dr = c1 - c2
        d2r = dr.sqr()
        R2ee += d2r 

        r_av = Real3D(0.,0.,0.)
        for k in range(mon_per_chain):
            pi = system.storage.getParticle(i*mon_per_chain+k+1)
            r_av += system.bc.getUnfoldedPosition(pi.pos, pi.imageBox)
        r_av /= mon_per_chain

        for k in range(mon_per_chain):
            pi = system.storage.getParticle(i*mon_per_chain+k+1)
            unf_i = system.bc.getUnfoldedPosition(pi.pos, pi.imageBox)
            dr = Real3D(unf_i.x - r_av.x, unf_i.y - r_av.y, unf_i.z - r_av.z)
            d2r = dr.sqr()
            R2gyr += d2r

    R2ee /= num_chains
    R2gyr /= (num_chains * mon_per_chain)
    s = str(R2ee)
    s2 = str(R2gyr)
    f_stat.write(s+'\t'+s2+'\n')
    f_stat.close()

main()
