#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

###########################################################################
#                                                                         #
#  This example creates for a potential file with tabulated entries     #
#                                                                         #
#  Note: also very nice for plotting potentials                           #
#                                                                         #
###########################################################################

import espresso

from espresso import Real3D

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


potLJ  = espresso.interaction.LennardJones(epsilon=1.0, sigma=1.0, shift=0.0, cutoff=2.5)
potFENE = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
potHarmonic = espresso.interaction.Harmonic(K=200.0, r0=1.2)
potCosine = espresso.interaction.Cosine(K=1.5, theta0=3.1415926)
potAngularHarmonic = espresso.interaction.AngularHarmonic(K=1.5, theta0=3.1415926)
potAngularCosineSquared = espresso.interaction.AngularCosineSquared(K=1.5, theta0=3.1415926)
potOPLS = espresso.interaction.OPLS(K1=1.0, K2=1.0, K3=1.0, K4=1.0)


#  Do not start at zero 
# 2-body
writeTabFile(potLJ, "pot-lj.txt", N=257, low=0.01, high=potLJ.cutoff)
writeTabFile(potFENE, "pot-fene.txt", N=257, low=0.0001, high=1.49)
writeTabFile(potHarmonic, "pot-harmonic.txt", N=512, low=0.001, high=15)

# 3-body
#writeTabFile(potCosine, "pot-cosine.txt", N=257, low=0.0001, high=3.14, body=3)
#writeTabFile(potAngularHarmonic, "pot-ah.txt", N=257, low=0.0001, high=3.14, body=3)
#writeTabFile(potAngularCosineSquared, "pot-acs.txt", N=257, low=0.0001, high=3.14, body=3)

# 4-body
#writeTabFile(potOPLS, "pot-opls.txt", N=128, low=-3.14, high=3.14, body=4)

