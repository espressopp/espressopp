import espresso
from espresso import Int3D, Real3D
from espresso.tools import decomp
from espresso.tools.init_cfg import lattice
import random
         
# some global definitions
rho          = 0.8442
num_particles= 20**3
skin         = 0.3
rc           = 2.5
epsilon      = 1.0
sigma        = 1.0
shift        = False
dt           = 0.01
gamma        = 1.0
temperature  = 1.0

ptrng=random
ptrng.seed(335977)

# Parallel Tempering (replica exchange) integrator
ptthermostats=[] 
pt = espresso.ParallelTempering(NumberOfSystems = 4, RNG = ptrng)
for i in range(0, pt.getNumberOfSystems()):
    pt.startDefiningSystem(i)
    
    random.seed(3456+i*43785)
    rand=random.random()
    x, y, z, Lx, Ly, Lz  = lattice.create(num_particles, rho, perfect=False, RNG=rand)
    boxsize              = (Lx, Ly, Lz)
    system               = espresso.System()
    rng                  = espresso.esutil.RNG()
    bc                   = espresso.bc.OrthorhombicBC(rng, boxsize)
    system.bc            = bc
    system.rng           = rng
    system.skin          = skin
    nodeGrid             = espresso.tools.decomp.nodeGrid(pt.getNumberOfCPUsPerSystem())
    cellGrid             = espresso.tools.decomp.cellGrid(boxsize,nodeGrid,rc,skin)
    storage              = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
    system.storage       = storage
    vl                   = espresso.VerletList(system,cutoff=rc+skin)
    potLJ                = espresso.interaction.LennardJones(epsilon, sigma, rc, shift)
    interLJ              = espresso.interaction.VerletListLennardJones(vl)
    integrator           = espresso.integrator.VelocityVerlet(system)
    integrator.dt        = dt
    langevin             = espresso.integrator.Langevin(system)
    langevin.gamma       = gamma
    langevin.temperature = temperature*i/10 + 0.5
    integrator.langevin  = langevin
    interLJ.setPotential(type1=0, type2=0, potential=potLJ)
    system.addInteraction(interLJ)
    for pid in range(num_particles):
        storage.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
    storage.decompose()

    pt.setIntegrator(integrator, langevin)
    pt.setAnalysisE(interLJ)
    pt.setAnalysisT(espresso.analysis.Temperature(system))
    pt.endDefiningSystem(i)

for p in range(10):
    pt.run(200)
    pt.exchange()
    multiT     = pt._multisystem.runAnalysisTemperature()
    print "%s" % multiT
