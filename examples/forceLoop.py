import espresso
import espresso.bc
import espresso.interaction
import espresso.pairs

# implied
import espresso.system

pbc=espresso.bc.PBC(length=10)
particles=espresso.system.Particles()

for pid in range(100):
    # default: pos=0 0 0
    particles.addParticle(pos=mySystem.randomPos())

allpairs=espresso.pairs.AllPairs(bc=pbc, group=particles)

# default: shift=auto, offset=0.0
ljint=espresso.interaction.LennardJones(sigma=1, epsilon=1, cutoff=2.5)

# vvintegrator=espresso.integrator.VelocityVerletIntegrator(timestep=0.001)
# vvintegrator.addForce(interaction=ljint, pairs=allpairs)
# vvintegrator.addForce(interaction=feneint, pairs=bonds)

allpairs.computeForces(setforce='force', interaction=ljint)

force=particles.getParticleProperty(name='force')

print force[particles[17]]

print particles[17].get(name='force')

# system stores global default variables 
#  * geometry
#  * pairs
#  * particles
