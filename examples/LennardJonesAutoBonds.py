import MPI
import logging
import espresso
from espresso import Real3D

box                  = (8,8,8)
rc                   = pow(2, 1.0/6.0)
rc_max               = 2.0
skin                 = 0.3
sigma                = 1.0
epsilon              = 1.0
nodeGrid             = espresso.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
cellGrid             = espresso.tools.decomp.cellGrid(box, nodeGrid, rc_max, skin)
system               = espresso.System()
system.rng           = espresso.esutil.RNG()
system.bc            = espresso.bc.OrthorhombicBC(system.rng, box)
system.skin          = skin
system.storage       = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

WCA_verlet_list      = espresso.VerletList(system, cutoff = rc + system.skin)
Harmonic_pair_bonds  = espresso.FixedPairList(system.storage)

system.storage.addParticles( [ [0, Real3D(4,1,1), Real3D(0,0,0)], \
                               [1, Real3D(5.4,1,1), Real3D(-0.3,0.3,0.0)] ], \
                             'id','pos','v')
system.storage.decompose()

WCA_potential        = espresso.interaction.LennardJonesAutoBonds(sigma, epsilon, cutoff = rc, bondlist = Harmonic_pair_bonds)
WCA_interaction      = espresso.interaction.VerletListLennardJonesAutoBonds(WCA_verlet_list)
WCA_interaction.setPotential(type1 = 0, type2 = 0, potential = WCA_potential)
Harmonic_potential   = espresso.interaction.Harmonic(K=10.0, r0=1.0)
Harmonic_interaction = espresso.interaction.FixedPairListHarmonic(system, Harmonic_pair_bonds, Harmonic_potential)

system.addInteraction(WCA_interaction)
system.addInteraction(Harmonic_interaction)

integrator           = espresso.integrator.VelocityVerlet(system)
integrator.dt        = 0.0001
sock = espresso.tools.vmd.connect(system)
# logging.getLogger("MDIntegrator").setLevel(logging.INFO) 
for i in range (1000000):
  integrator.run(100)
  # print 'position=',system.storage.getParticle(0).pos, system.storage.getParticle(1).pos
  # print Harmonic_pair_bonds.getBonds()
  espresso.tools.vmd.imd_positions(system, sock)

# logging.getLogger("MDIntegrator").setLevel(logging.FATAL) 


