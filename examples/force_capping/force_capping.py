import espressopp

# try to run this script with force_capping switched on and off
forceCapping = True

# specify number of particles
num_particles = 200
rho = 0.85
L   = pow(num_particles/rho, 1.0/3.0)

# setup random LennardJones system
# this system is likely to explode on integration, because particles can strongly overlap
system, integrator = espressopp.standard_system.LennardJones(num_particles, box=(L, L, L), temperature=1.0)

# choose a smaller timestep
integrator.dt = 0.0001

if forceCapping:
  max_force = 100000.0
  # define force capping extension
  capForce = espressopp.integrator.CapForce(system, max_force)
  # and add it to the integrator
  integrator.addExtension(capForce)

espressopp.tools.analyse.info(system, integrator)

sock = espressopp.tools.vmd.connect(system)
for i in range(1000):
  # make 10 Velocity-Verlet integration steps
  integrator.run(10)
  # print system information
  espressopp.tools.analyse.info(system, integrator)
  # update postions in VMD
  espressopp.tools.vmd.imd_positions(system, sock)
  #switch off force capping after some time
  if forceCapping and i>100:
    capForce.disconnect()
    forceCapping = False
    print "switching off force capping"

