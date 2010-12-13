import espresso
from espresso import Int3D, Real3D
from espresso.tools.init_cfg import lattice

rho          = 0.8442
num_particles= 20**3
skin         = 0.3
rc           = 2.5
epsilon      = 1.0
sigma        = 1.0
shift        = False
dt           = 0.005
gamma        = 10.0
temperature1 = 1.0
temperature2 = 0.5
nodeGrid     = Int3D(2,1,1)
cellGrid     = Int3D(3,3,3)
rng=espresso.esutil.RNG()
x, y, z, Lx, Ly, Lz = lattice.create(num_particles, rho, perfect=False, RNG=rng)
size = (Lx, Ly, Lz)

multisystem = espresso.MultiSystem()

################################################
# Setup system1
################################################
comm1=espresso.pmi.Communicator([0,1])
espresso.pmi.activate(comm1)
multisystem.beginSystemDefinition()

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
for pid in range(num_particles):
    storage1.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
storage1.decompose()
integrator1           = espresso.integrator.VelocityVerlet(system1)
integrator1.dt        = dt
langevin1             = espresso.integrator.Langevin(system1)
langevin1.gamma       = gamma
langevin1.temperature = temperature1
integrator1.langevin  = langevin1
analysisT1 = espresso.analysis.Temperature(system1)

multisystem.setIntegrator(integrator1)
multisystem.setAnalysisPotential(interLJ1)
multisystem.setAnalysisTemperature(analysisT1)
espresso.pmi.deactivate(comm1)

################################################
# Setup system2
################################################
comm2=espresso.pmi.Communicator([2,3])
espresso.pmi.activate(comm2)
multisystem.beginSystemDefinition()

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
for pid in range(num_particles):
    storage2.addParticle(pid + 1, Real3D(x[pid], y[pid], z[pid]))
storage2.decompose()
integrator2           = espresso.integrator.VelocityVerlet(system2)
integrator2.dt        = dt
langevin2             = espresso.integrator.Langevin(system2)
langevin2.gamma       = gamma
langevin2.temperature = temperature2
integrator2.langevin  = langevin2
analysisT2 = espresso.analysis.Temperature(system2)

multisystem.setIntegrator(integrator2)
multisystem.setAnalysisPotential(interLJ2)
multisystem.setAnalysisTemperature(analysisT2)
espresso.pmi.deactivate(comm2)

print "Potential Energy of system1 is ", multisystem.runAnalysisPotential()[0]
print "Potential Energy of system2 is ", multisystem.runAnalysisPotential()[2]
print "Temperature of system1 is ", multisystem.runAnalysisTemperature()[0]
print "Temperature of system2 is ", multisystem.runAnalysisTemperature()[2]

multisystem.runIntegrator(100)

print "Potential Energy of system1 is ", multisystem.runAnalysisPotential()[0]
print "Potential Energy of system2 is ", multisystem.runAnalysisPotential()[2]
print "Temperature of system1 is ", multisystem.runAnalysisTemperature()[0]
print "Temperature of system2 is ", multisystem.runAnalysisTemperature()[2]
