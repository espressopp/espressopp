import espresso
from espresso import Int3D, Real3D
from espresso.tools import decomp
from espresso.tools.init_cfg import lattice
import random
import sys
         
# some global definitions
skin         = 0.3
rc           = 2.5
epsilon      = 1.0
sigma        = 1.0
shift        = 0.0
dt           = 0.005
gamma        = 1.0
temperature  = 1.0

ptrng=random
ptrng.seed(335977)

if espresso.MPI.COMM_WORLD.size != 4:
  print "currently this example can only be run with 4 CPUs"
  sys.exit(0)

# Parallel Tempering (replica exchange) integrator
ptthermostats=[] 
pt = espresso.ParallelTempering(NumberOfSystems = 4, RNG = ptrng)
for i in range(0, pt.getNumberOfSystems()):
    pt.startDefiningSystem(i)
    pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = espresso.tools.readxyz('parallel_tempering.xyz')
    num_particles        = len(pid)
    boxsize              = (Lx, Ly, Lz)
    rho                  = num_particles / (Lx * Ly * Lz)
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
    vl                   = espresso.VerletList(system,cutoff=rc)
    potLJ                = espresso.interaction.LennardJones(epsilon, sigma, rc, shift)
    interLJ              = espresso.interaction.VerletListLennardJones(vl)
    integrator           = espresso.integrator.VelocityVerlet(system)
    integrator.dt        = dt
    langevin             = espresso.integrator.LangevinThermostat(system)
    langevin.gamma       = gamma
    langevin.temperature = temperature*i/10 + 0.5
    integrator.addExtension(langevin)
    interLJ.setPotential(type1=0, type2=0, potential=potLJ)
    system.addInteraction(interLJ)
    for k in range(num_particles):
        storage.addParticle(pid[k], Real3D(x[k], y[k], z[k]), checkexist=False)
    storage.decompose()

    pt.setIntegrator(integrator, langevin)
    pt.setAnalysisE(interLJ)
    pt.setAnalysisT(espresso.analysis.Temperature(system))
    pt.setAnalysisNPart(espresso.analysis.NPart(system))
    pt.endDefiningSystem(i)

for p in range(10):
    pt.run(200)
    pt.exchange()
    multiT     = pt._multisystem.runAnalysisTemperature()
    print "%s" % multiT
