import espresso
import MPI

def Default(box, rc=1.12246, skin=0.3, dt=0.005, temperature=None):
  '''
  return default system and integrator, no interactions, no particles are set
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
    thermostat             = espresso.integrator.Langevin(system)
    thermostat.gamma       = 1.0
    thermostat.temperature = temperature
    integrator.langevin    = thermostat
   
  return system, integrator
