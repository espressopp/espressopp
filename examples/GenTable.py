#!/usr/bin/env python

###########################################################################
#                                                                         #
#  This example creates for a potential a file with tabulated entries     #
#                                                                         #
#  Note: also very nice for plotting potentials                           #
#                                                                         #
###########################################################################

import espresso

from espresso import Real3D

def writeTabFile(pot, name, N, low = 0.0):

    outfile = open(name, "w")

    high = pot.cutoff

    delta = (high - low) / (N - 1)

    for i in range(N):

        r = low + i * delta;
        energy = pot.computeEnergy(r)

        dist = Real3D(r, 0.0, 0.0)
        force = pot.computeForce(dist)
        force = force[0]

        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))

    outfile.close()

cutoff = 3.0

potLJ  = espresso.interaction.LennardJones(sigma = 1.0, epsilon = 1.0, shift = 0.0, cutoff = cutoff)

#  Do not start at zero 

writeTabFile(potLJ, "pair.txt", N = 257, low = 0.5)
