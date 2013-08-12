import espresso

system, integrator = espresso.standard_system.LennardJones(0,(10,10,10))
system.storage.addParticles([[1,0,espresso.Real3D(5,5,5)]],'id','type','pos')
system.storage.addParticles([[2,0,espresso.Real3D(5,6,5)]],'id','type','pos')
system.storage.addParticles([[3,1,espresso.Real3D(6,5,5)]],'id','type','pos')
system.storage.addParticles([[4,1,espresso.Real3D(6,6,5)]],'id','type','pos')
system.storage.decompose()

Epot = espresso.analysis.EnergyPot(system)
print 'Epot = ',Epot.compute()

intLJ = system.getInteraction(0)
pot00 = intLJ.clonePotential(0,0)
intLJ.setPotential(1,1,pot00)
intLJ.setPotential(1,0,pot00)
intLJ.setPotential(0,1,pot00)

print 'Epot = ', Epot.compute()

