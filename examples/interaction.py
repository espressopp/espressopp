import espresso.interaction

lj = espresso.interaction.LennardJones(epsilon=2.0)
print("lj(1.2) = %f" % lj.computeEnergy(1.2))
print("lj(1.3) = %f" % lj.computeEnergy(1.3))

fene = espresso.interaction.FENE(K=1.5, r0=1.0, rMax=1.0)
print("fene(1.2) = %f" % fene.computeEnergy(1.2))
print("fene(1.3) = %f" % fene.computeEnergy(1.3))



