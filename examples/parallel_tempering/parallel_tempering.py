import espressopp
from espressopp import Int3D, Real3D
from espressopp.tools import decomp
from espressopp.tools import lattice
import random
import sys

print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print "+ Due to a recent change in the design of ESPResSo++ the Parallel Tempering     +"
print "+ Class is not available in the current version. This Class will be back soon.  +"
print "+ Multisystem simulations are still possible but have to be setup manually.     +"
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

sys.exit(0) 
         
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

if espressopp.MPI.COMM_WORLD.size != 4:
  print "currently this example can only be run with 4 CPUs"
  sys.exit(0)

# Parallel Tempering (replica exchange) integrator
ptthermostats=[] 
pt = espressopp.ParallelTempering(NumberOfSystems = 4, RNG = ptrng)
for i in range(0, pt.getNumberOfSystems()):
    pt.startDefiningSystem(i)
    pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz('parallel_tempering.xyz')
    num_particles        = len(pid)
    boxsize              = (Lx, Ly, Lz)
    rho                  = num_particles / (Lx * Ly * Lz)
    system               = espressopp.System()
    rng                  = espressopp.esutil.RNG()
    bc                   = espressopp.bc.OrthorhombicBC(rng, boxsize)
    system.bc            = bc
    system.rng           = rng
    system.skin          = skin
    nodeGrid             = espressopp.tools.decomp.nodeGrid(pt.getNumberOfCPUsPerSystem())
    cellGrid             = espressopp.tools.decomp.cellGrid(boxsize,nodeGrid,rc,skin)
    storage              = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid, nocheck=True)
    system.storage       = storage
    vl                   = espressopp.VerletList(system,cutoff=rc)
    potLJ                = espressopp.interaction.LennardJones(epsilon, sigma, rc, shift)
    interLJ              = espressopp.interaction.VerletListLennardJones(vl)
    integrator           = espressopp.integrator.VelocityVerlet(system)
    integrator.dt        = dt
    langevin             = espressopp.integrator.LangevinThermostat(system)
    langevin.gamma       = gamma
    langevin.temperature = temperature*i/20 + 0.2
    integrator.addExtension(langevin)
    interLJ.setPotential(type1=0, type2=0, potential=potLJ)
    system.addInteraction(interLJ)
    for k in range(num_particles):
        storage.addParticle(pid[k], Real3D(x[k], y[k], z[k]), checkexist=False)
    storage.decompose()

    pt.setIntegrator(integrator, langevin)
    pt.setAnalysisE(interLJ)
    pt.setAnalysisT(espressopp.analysis.Temperature(system))
    pt.setAnalysisNPart(espressopp.analysis.NPart(system))
    pt.endDefiningSystem(i)

# let each system reach its temperature
for p in range(100):
    pt.run(100)
    multiT     = pt._multisystem.runAnalysisTemperature()
    print "%s" % multiT

for p in range(10):
    pt.run(200)
    pt.exchange()
    multiT     = pt._multisystem.runAnalysisTemperature()
    print "%s" % multiT
