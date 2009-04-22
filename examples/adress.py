import espresso

# Parameters of ADRes
# * weighting function
# * function to add details of the different resolution models
# * 

# * we will have a list of pairs of a particle in one resolution model
#   and a particle in the other resolution model

# An adress scheme needs to provide:
# * A function that can be hooked to
#   * the force computations (to exchange forces)
#   * the position update (to adapt the position of the CG model)
#   * at the end of the intergration (to add/remove details in one or
#     the other model)

# Integration:
# * computes forces
# * updates positions
# * communicate positions

pbc = espresso.bc.PBC(length=10)

# SET UP THE COARSE-GRAINED REGION
cgParticles = espresso.decomposition.CellStorage(bc=pbc, grid=(2,2,2), skin=0.1)
cgIntegrator = espresso.integrator.VelocityVerlet(timestep=0.005, pb=pbc)
cgLJInt = espresso.interaction.LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.0)
cgVerletLists = espresso.pairs.VerletLists(radius=cgLJInt.cutoff, skin=0.3)
cgIntegrator.addInteraction(pairs=cgVerletLists, interaction=cgLJInt)

# SET UP THE ATOMIC SYSTEM
atomicParticles = espresso.decomposition.CellStorage(bc=pbc, grid=(2,2,2), skin=0.1)
atomicIntegrator = espresso.integrator.VelocityVerlet(timestep=0.001, bc=pbc)
atomicBonds = espresso.pairs.List()
atomicFENEInt = espresso.interaction.FENE(r0=0.5, K=1.0, rMax=0.1)
atomicIntegrator.addInteraction(pairs=atomicBonds, interaction=atomicFENEInt)

# Python class that creates the atomic resolution particles for a
# given coarse-grained particle
class MyCreator(object):
    def __init__(self, atomicParticles, atomicBonds):
        # store references to required objects
        self.cgParticles=cgParticles
        self.atomicBonds=atomicBonds
        # init an RNG
        import random
        self.rnd = random.Random()
        # set up a library of example hires particles
        self.createLib()

    def create(self, cgParticle):
        # get the position of the particle
        pos=cgParticle[pos]
        # draw a random example from the library
        pos1, pos2 = self.rnd.sample(self.library)
        # add the new atomic particles 
        p1=self.atomicParticles.addParticle(pos+pos1)
        p2=self.atomicParticles.addParticle(pos+pos2)
        # create the bond between them
        self.atomicBonds.add(p1, p2)
        # and return the particles that were created
        return p1, p2

# create an instance of the creator
myCreator = MyCreator(atomicParticles, atomicBonds)

# DEFINE THE ADRESS INTEGRATOR
adInt = espresso.adress.Integrator(
    lowresIntegrator=cgIntegrator,
    lowresParticles=cgParticles,
    hiresIntegrator=atomicIntegrator,
    hiresParticles=atomicParticles,
    weighting=espresso.adress.SlabWeighting(
        normal='x',
        position=(5.0, 5.0, 5.0),
        size=2.0,
        hybridSize=1.0,
        createSkinSize=0.5,
        createTriggerSize=0.1,
        destroySkinSize=0.5,
        destroyTriggerSize=0.1),
    coupling=espresso.adress.StandardCoupling(),
    hiresCreator=myCreator)

# CREATE THE INITIAL PARTICLES
# .
# .

# SIMULATE
adInt.integrate(steps=1000)

#
# C++ weighting interface:
# * real getWeight(ParticleId &p)
# * bool triggersDestruction(ParticleId &p)
#   * returns whether the particle has atomic data and is in the destruction trigger skin
# * bool needsDestruction(ParticleId &p)
#   * returns whether the particle has atomic data and is in the destruction skin
# * bool triggersCreation(ParticleId &p)
#   * returns whether the particle has no atomic data and is in the construction trigger skin
# * bool needsCreation(ParticleId &p)
#   * returns whether the particle has no atomic data and is in the construction skin
# * adress.SlabWeighting:
#   * defines a slab region (part of a triclinic box)
#

#
# C++ coupling interface:
# * template <class ParticleIdContainer> coupleForces(ParticleId, ParticleIdContainer&)
# * template <class ParticleIdContainer> couplePositions(ParticleId, ParticleIdContainer&)
# * StandardCoupling:
#   * coupleForces
#     * will reweight the forces on the CG particle
#     * will reweight the forces on the AA particles
#     * will add the forces of the AA particles to the CG particles
#     * will distribute the force of the CG particle onto the AA particles
#   * couplePositions
#     * will set the position of the CG particle to the COM of the AA particles
#

#
# hiresCreator Python interface:
# hiresCreator.((Particle)p) --> [ Particle ]
#   * creates the high resolution information for the cg particle p1, returns the list of created particles
# 

#
# adress.Integrator: each step, it will be checked whether a particle has moved further into the skin
