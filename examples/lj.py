import espresso._espresso

lj = espresso._espresso.interaction_LennardJones()

print "lj(1.2) =", lj.computeEnergy(1.2)
print "lj(1.3) =", lj.computeEnergy(1.3)
