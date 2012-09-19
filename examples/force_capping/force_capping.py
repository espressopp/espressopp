import espresso

# try to run this script with force_capping switched on and off
forceCapping = True

# specify number of particles
num_particles = 200
rho = 0.85
L   = pow(num_particles/rho, 1.0/3.0)

# setup random LennardJones system
# this system is likely to explode on integration, because particles can strongly overlap
system, integrator = espresso.standard_system.LennardJones(num_particles, box=(L, L, L), temperature=1.0)

# use a very small timestep
integrator.dt = 0.0001

if forceCapping:
  max_force = 100000.0
  # define force capping extension
  capForce = espresso.integrator.CapForce(system, max_force)
  # and add it to the integrator
  integrator.addExtension(capForce)

espresso.tools.analyse.info(system, integrator)

integrator.run(50)
capForce.disconnect()

sock = espresso.tools.vmd.connect(system)
for i in range(1000):
  # make 10 Velocity-Verlet integration steps
  integrator.run(100)
  espresso.tools.analyse.info(system, integrator)
  #switch off force capping after some time
  if forceCapping and i>1:
    capForce.disconnect()
    forceCapping = False
  # update postions in VMD
  espresso.tools.vmd.imd_positions(system, sock)

