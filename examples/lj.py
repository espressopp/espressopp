# import espresso._espresso

import espresso._espresso

lj = espresso._espresso.interaction_LennardJones()

print "lj(1.2) =", lj.computeEnergy(1.2)
print "lj(1.3) =", lj.computeEnergy(1.3)

fene = espresso._espresso.interaction_FENE()

fene.setK(1.5)
fene.setr0(1.0)
fene.setrMax(1.0)

print "fene(1.2) =", fene.computeEnergy(1.2)
print "fene(1.3) =", fene.computeEnergy(1.3)

