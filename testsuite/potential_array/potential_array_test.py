#!/usr/bin/env python2

import espressopp
from espressopp.tools import replicate, readxyz
from espressopp import Real3D

import time, sys
import numpy as np

def bench_lj(xyz_file):

    # read from xyz file
    print "reading from: ", xyz_file
    pid, type, xpos, ypos, zpos, xvel, yvel, zvel, Lx, Ly, Lz = readxyz(xyz_file)
    Npart              = len(pid)
    box                = (Lx,Ly,Lz)

    ################################################
    # SIMULATION PARAMETERS
    ################################################

    # cutoff of the short range potential
    r_cutoff           = 2.5
    # VerletList skin size (also used for domain decomposition)
    skin               = 0.4
    # time step for the velocity verlet integrator
    dt                 = 0.005
    # Lennard Jones epsilon during equilibration phase
    epsilon00          = 0.1
    epsilon11          = 0.9
    epsilon01          = (0.1*0.3)**(0.5)
    # Lennard Jones sigma during warmup and equilibration
    sigma00            = 1.0
    sigma11            = 0.6
    sigma01            = (sigma00+sigma11)/2.0

    steps              = 100

    ################################################
    # SETUP SYSTEM
    ################################################

    # create the basic system
    system             = espressopp.System()
    # use the random number generator that is included within the ESPResSo++ package
    system.rng         = espressopp.esutil.RNG()
    # use orthorhombic periodic boundary conditions
    system.bc          = espressopp.bc.OrthorhombicBC(system.rng, box)
    # set the skin size used for verlet lists and cell sizes
    system.skin        = skin
    # get the number of CPUs to use
    NCPUs              = espressopp.MPI.COMM_WORLD.size
    # calculate a regular 3D grid according to the number of CPUs available
    nodeGrid           = espressopp.tools.decomp.nodeGrid(NCPUs, box, r_cutoff, skin)
    # calculate a 3D subgrid to speed up verlet list builds and communication
    cellGrid           = espressopp.tools.decomp.cellGrid(box, nodeGrid,  r_cutoff, skin)

    print "gitrevision        = ", espressopp.Version().gitrevision
    print "NCPUs              = ", NCPUs
    print "r_cutoff           = ", r_cutoff
    print "skin               = ", skin
    print "box                = ", box
    print "nodeGrid           = ", nodeGrid
    print "cellGrid           = ", cellGrid
    print "Nparticles         = ", Npart
    print "steps              = ", steps

    # create a domain decomposition particle storage with the calculated nodeGrid and cellGrid
    system.storage     = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    # use a velocity Verlet integration scheme
    integrator         = espressopp.integrator.VelocityVerlet(system)

    # set the integration step
    integrator.dt  = dt


    ################################################
    # ADD PARTICLES
    ################################################

    props = ['id', 'type', 'pos']

    print("adding particles to storage...")
    new_particles = []
    pid = 0
    for x, y, z in zip(xpos, ypos, zpos):
        ptype = pid % 2
        new_particles.append([pid, ptype, Real3D(x,y,z)])
        pid += 1
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()

    print("adding particles to storage... DONE")

    ########################################################################
    # 7. setting up interaction potential for the equilibration            #
    ########################################################################
    verletlist  = espressopp.VerletList(system, r_cutoff)
    interaction = espressopp.interaction.VerletListLennardJones(verletlist)
    interaction.setPotential(type1=0, type2=0,
                    potential=espressopp.interaction.LennardJones(
                    epsilon=epsilon00, sigma=sigma00, cutoff=r_cutoff, shift=0.0))
    interaction.setPotential(type1=0, type2=1,
                    potential=espressopp.interaction.LennardJones(
                    epsilon=epsilon01, sigma=sigma01, cutoff=r_cutoff, shift=0.0))
    interaction.setPotential(type1=1, type2=1,
                    potential=espressopp.interaction.LennardJones(
                    epsilon=epsilon11, sigma=sigma11, cutoff=r_cutoff, shift=0.0))
    system.addInteraction(interaction)

    print("running integrator...")
    integrator.run(steps)
    print("running integrator... DONE")

    keys = [
        "Run",                  # 0
        "ForceComp[0]",         # 1
        "ForceComp[1]",         # 2
        "ForceComp[2]",         # 3
        "UpdateGhosts",         # 4
        "CollectGhostForces",   # 5
        "Integrate1",           # 6
        "Integrate2",           # 7
        "Resort",               # 8
        "Lost"                  # 9
    ]
    sub_keys = [0,6,8,4,1,2,5,7]

    timers = list(integrator.getTimers())
    timers_summ = np.mean(np.array(timers),axis=0)

    # add non-interacting particles of different type
    # this forces potentialArray to resize
    print("adding more particles to storage...")
    new_particles = []
    for x, y, z in zip(xpos[:10], ypos[:10], zpos[:10]):
        ptype = 6
        new_particles.append([pid, ptype, Real3D(x,y,z)])
        pid += 1
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()

    print("adding more particles to storage... DONE")

    print("running integrator...")
    integrator.run(steps)
    print("running integrator... DONE")

    timers = list(integrator.getTimers())
    timers_summ += np.mean(np.array(timers),axis=0)

    print ""
    for k in sub_keys:
        print "{:18} =  {}".format(keys[k],timers_summ[k])
    print ""

    temperature = espressopp.analysis.Temperature(system)
    T = temperature.compute()
    print "temperature        =  {}".format(T)

    T_exp = 0.160074388303
    assert(abs(T-T_exp)<1.0e-8)

def main():
    xyz_file = "lennard_jones_fluid_10000_2048.xyz"
    bench_lj(xyz_file)

if __name__ == '__main__':
    main()
