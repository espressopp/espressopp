import espresso
from espresso import Int3D

L            = 10
density      = 0.01
size         = (L, L, L)
npart        = int(L * L * L * density)
skin         = 0.3
rc           = 2.5
epsilon      = 1.0
sigma        = 1.0
shift        = False
dt           = 0.000000001
gamma        = 10.0
temperature1 = 1.0
temperature2 = 0.5
nodeGrid     = Int3D(2,1,1)
cellGrid     = Int3D(3,3,3)

def create_VV_integrator_with_Langevin(s,dt,gamma,temperature) :
    integrator = espresso.integrator.VelocityVerlet(s)
    integrator.dt = dt
    integrator.langevin = espresso.integrator.Langevin(s)
    integrator.langevin.gamma = gamma
    integrator.langevin.temperature = temperature
    return integrator


################################################
# Setup system1
################################################
comm1=espresso.pmi.Communicator([1,3])
espresso.pmi.activate(comm1)

system1         = espresso.System()
rng1            = espresso.esutil.RNG()
bc1             = espresso.bc.OrthorhombicBC(rng1, size)
system1.bc      = bc1
system1.rng     = rng1
system1.skin    = skin
storage1        = espresso.storage.DomainDecomposition(system1, nodeGrid, cellGrid)
system1.storage = storage1
vl1             = espresso.VerletList(system1,cutoff=rc+skin)
potLJ1          = espresso.interaction.LennardJones(epsilon, sigma, rc, shift)
interLJ1        = espresso.interaction.VerletListLennardJones(vl1)
interLJ1.setPotential(type1=0, type2=0, potential=potLJ1)
system1.addInteraction(interLJ1)
for pid in range(npart) :
    print "system1:",pid
    storage1.addParticle(pid+1, bc1.getRandomPos())
storage1.decompose()

espresso.pmi.deactivate(comm1)

################################################
# Setup system2
################################################
comm2=espresso.pmi.Communicator([0,2])
espresso.pmi.activate(comm2)

system2         = espresso.System()
rng2            = espresso.esutil.RNG()
bc2             = espresso.bc.OrthorhombicBC(rng2, size)
system2.bc      = bc2
system2.rng     = rng2
system2.skin    = skin
storage2        = espresso.storage.DomainDecomposition(system2, nodeGrid, cellGrid)
system2.storage = storage2
vl2             = espresso.VerletList(system2,cutoff=rc+skin)
potLJ2          = espresso.interaction.LennardJones(epsilon, sigma, rc, shift)
interLJ2        = espresso.interaction.VerletListLennardJones(vl2)
interLJ2.setPotential(type1=0, type2=0, potential=potLJ2)
system2.addInteraction(interLJ2)
for pid in range(npart) :
    print "system2:",pid
    storage2.addParticle(pid+1, bc2.getRandomPos())
#storage2.decompose()

espresso.pmi.deactivate(comm2)


