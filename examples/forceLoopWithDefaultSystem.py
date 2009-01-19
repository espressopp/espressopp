import espresso
import espresso.bc
import espresso.interaction
import espresso.pairs

# implied
import espresso.system

pbc=espresso.bc.PBC(length=10)

mySystem=espresso.system.default

mySystem.set(bc=pbc)

for pid in range(100):
    # default: pos=0 0 0
    mySystem.addParticle(pos=mySystem.randomPos())

# implies that all particles are used?
# default: system=espresso.system.default, bc=$system.bc, group=$system.particles
allpairs=espresso.pairs.AllPairs()

# default: shift=auto, offset=0.0
ljint=espresso.interaction.LennardJones(sigma=1, epsilon=1, cutoff=2.5)

# vvintegrator=espresso.integrator.VelocityVerletIntegrator(timestep=0.001)
# vvintegrator.addForce(interaction=ljint, pairs=allpairs)
# vvintegrator.addForce(interaction=feneint, pairs=bonds)
# mySystem.set(integrator=vvintegrator)

ljint.computeForces(setforce='force',pairs=allpairs)

force=mySystem.getParticleProperty(name='force')
print force[mySystem.particles[17]]

print mySystem.particles[17].get(name='force')

# system stores global default variables 
#  * geometry
#  * pairs
#  * particles
