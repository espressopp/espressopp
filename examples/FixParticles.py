import espresso

rc             = pow(2.0, 1.0/6.0)
box            = (10, 10, 10)
dt             = 0.005
skin           = 0.3
epsilon        = 1.0
sigma          = 1.0
temperature    = None
system         = espresso.System()
system.rng     = espresso.esutil.RNG()
system.bc      = espresso.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
nodeGrid       = espresso.tools.decomp.nodeGrid(espresso.MPI.COMM_WORLD.size)
cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
interaction    = espresso.interaction.VerletListLennardJones(espresso.VerletList(system, cutoff=rc))
interaction.setPotential(type1=0, type2=0, potential=espresso.interaction.LennardJones(epsilon, sigma, rc))
system.addInteraction(interaction)

c_fixed = 0
c_free  = 1
# fix all 3 coordinates
fixMask = espresso.Int3D(c_fixed, c_fixed, c_fixed)

# create a particel group that will contain the fixed particles
fixedWall  = espresso.ParticleGroup(system.storage)

# add a particle wall
pid = 1
for k in range(10):
  for l in range(10):
    system.storage.addParticle(pid, espresso.Real3D(k, 5, l))
    fixedWall.add(pid)
    pid += 1

# add also one free particle
system.storage.addParticle(0, espresso.Real3D(5.8,9,5.5))
system.storage.modifyParticle(0, 'v', espresso.Real3D(0, -0.1, 0))

# don't forget do decompose !
system.storage.decompose()

# setup the integrator
# here we want to use the special VelocityVerletFixedParticle integrator
integrator     = espresso.integrator.VelocityVerletFixedParticles(system)
integrator.fixedParticles = fixedWall
integrator.setFixMask(fixMask)
integrator.dt  = dt
if (temperature != None):
  langevin             = espresso.integrator.LangevinThermostat(system)
  langevin.gamma       = 1.0
  langevin.temperature = temperature
  integrator.addExtension(langevin)

sock = espresso.tools.vmd.connect(system)
for i in range(10000):
  integrator.run(100)
  espresso.tools.vmd.imd_positions(system, sock)
