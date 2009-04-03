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

pbc=espresso.bc.PBC(length=10)

atomic=espresso.decomposition.AtomicStorage()
# define the rest of the atomic region: forces, etc.
atomicIntegrator = espresso.integrator.VelocityVerlet(timestep=0.001, pb=pbc)

coarse = espresso.decomposition.AtomicStorage()
# define the rest of the cg region: forces etc.
coarseIntegrator = espresso.integrator.VelocityVerlet(timestep=0.005, pb=pbc)

adInt = adress.Integrator(
    lowresIntegrator=coarseIntegrator,
    hiresIntegrator=atomicIntegrator,
    weighting=adress.SlabWeighting(normal='x',
        position=(5.0, 5.0, 5.0),
        size=2.0,
        hybridSize=1.0,
        createSkinSize=0.5,
        createTriggerSize=0.1,
        destroySkinSize=0.5,
        destroyTriggerSize=0.1),
    coupling=adress.StandardCoupling(),
    hiresCreator=myHiresCreator)

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

adInt.integrate(steps=1000)
