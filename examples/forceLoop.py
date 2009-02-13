import sys
import espresso
import espresso.bc
import espresso.interaction
import espresso.pairs

pbc=espresso.bc.PBC(length=10)
particles=espresso.decomposition.AtomicStorage()
#particles=espresso.decomposition.CellStorage(bc=pbc)

for pid in range(100):
    # default: pos=0 0 0
    particles.addParticle(pos=pbc.randomPos())

allpairs=espresso.pairs.All(bc=pbc,set=particles)

# default: shift=auto, offset=0.0
ljint=espresso.interaction.LennardJones(sigma=1, epsilon=1, cutoff=2.5)

# Just compute the forces
allpairs.computeForces(setforce='force', interaction=ljint)
force=particles.getParticleProperty(name='force')
sys.write(force[particles[17]])
sys.write(particles[17].get(name='force'))

# Now do 1000 simulation steps
vvintegrator=espresso.integrator.VelocityVerlet(timestep=0.001, bc=pbc)
vvintegrator.addForce(interaction=ljint, pairs=allpairs)
vvintegrator.run(steps=1000)

#vvintegrator.addForce(interaction=feneint, pairs=bonds)

# system stores global default variables 
#  * geometry
#  * pairs
#  * particles
