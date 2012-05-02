import espresso
system, integrator = espresso.standard_system.LennardJones(0, (10, 10, 10))

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

# set the movement constraint of the integrator
integrator.fixpositions = espresso.integrator.FixPositions(system, fixedWall, fixMask)
integrator.dt = 0.01

sock = espresso.tools.vmd.connect(system)
for i in range(10000):
  integrator.run(100)
  espresso.tools.vmd.imd_positions(system, sock)
