"""
************************************
**espresso.standard_system.Minimal**
************************************

"""
import espresso
import MPI

def Minimal(num_particles, box, rc=1.12246, skin=0.3, dt=0.005, temperature=None):
  '''return minimal system and integrator whithout any interactions defined:
  particles have random positions in box
  if tempearture is != None then Langevin thermostat is set to temperature (gamma is 1.0)
  '''
  system         = espresso.System()
  system.rng     = espresso.esutil.RNG()
  system.bc      = espresso.bc.OrthorhombicBC(system.rng, box)
  system.skin    = skin
  nodeGrid       = espresso.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
  cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
  system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)

  integrator     = espresso.integrator.VelocityVerlet(system)  
  integrator.dt  = dt
  if (temperature != None):
    thermostat             = espresso.integrator.LangevinThermostat(system)
    thermostat.gamma       = 1.0
    thermostat.temperature = temperature
    integrator.addExtension(thermostat)
  
  props = ['id', 'type', 'mass', 'pos', 'v']
  new_particles = []
  pid = 1
  while pid <= num_particles:
    type = 0
    mass = 1.0
    pos  = system.bc.getRandomPos()
    vel  = espresso.Real3D(0.0, 0.0, 0.0)
    part = [pid, type, mass, pos, vel]
    new_particles.append(part)
    if pid % 1000 == 0:
        system.storage.addParticles(new_particles, *props)
        system.storage.decompose()
        new_particles = []
    pid += 1
  system.storage.addParticles(new_particles, *props)
  system.storage.decompose()
 
  return system, integrator
