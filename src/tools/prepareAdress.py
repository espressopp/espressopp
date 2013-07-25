"""

********************
**prepareAdress.py** 
********************

Auxiliary python functions for preparation of an Adress Simulation based on a configuration from an all-atomistic simulation.

If one uses a configuration file from an all-atomistic simulation as start configuration for an AdResS simulation,
the particles are probably all located inside the simulation box. However, in AdResS only the coarse-grained center-of-mass
particles have to be in the box, the atomistic particles of the coarse grained might be outside around their CoM CG particle. When in the start configuration
atomistic particles belonging to a molecule are folded such that some of the atoms are on the one side of the box while the others are folded
to the other side the calculation of the center of mass goes wrong and the simulation will be incorrect. This script ensures a proper
center of mass calculation and a proper folding and configuration for the AdResS simulation by simply putting the CG particle in one of the atoms (AdressSetCG) first.
Then the molecules will be put together properly afterwards when calling AdressDecomp.

"""

import espresso

def AdressSetCG(apm, pidCG, allParticlesAT):
    cmp = [0,0,0]
    pos = (allParticlesAT[pidCG*apm])[1]
    for i in range(3):
	cmp[i] += pos[i]
    return cmp

def AdressDecomp(system, integrator):
    #espresso.system.storage.decompose()
    integrator.run(0)
    system.storage.decompose()
    integrator.run(0)