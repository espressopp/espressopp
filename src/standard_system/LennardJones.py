"""
*****************************************
**espresso.standard_system.LennardJones**
*****************************************

"""
import espresso
import MPI

def LennardJones(num_particles, box=(0,0,0), rc=1.12246, skin=0.3, dt=0.005, epsilon=1.0, sigma=1.0, shift='auto', temperature=None, xyzfilename=None, xyzrfilename=None):
  '''return random Lennard Jones system and integrator:
  if tempearture is != None then Langevin thermostat is set to temperature (gamma is 1.0)
  '''
    
  if xyzfilename and xyzrfilename:
     print "ERROR: only one of xyzfilename (only xyz data) or xyzrfilename (additional particle radius data) can be provided."
     sys.exit(1)

  if xyzrfilename: 
    pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf, radiusf = espresso.tools.readxyzr(xyzrfilename)
    box = (Lxf, Lyf, Lzf)
    num_particles = len(pidf)
  elif xyzfilename:
    pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf = espresso.tools.readxyz(xyzfilename)
    box = (Lxf, Lyf, Lzf)
    num_particles = len(pidf)
  else:
    if box[0]<=0 or box[1]<=0 or box[2]<=0:
      print "WARNING: no valid box size specified, box size set to (100,100,100) !"    
    
  system         = espresso.System()
  system.rng     = espresso.esutil.RNG()
  system.bc      = espresso.bc.OrthorhombicBC(system.rng, box)
  system.skin    = skin
  nodeGrid       = espresso.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
  cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
  system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
  interaction    = espresso.interaction.VerletListLennardJones(espresso.VerletList(system, cutoff=rc))
  interaction.setPotential(type1=0, type2=0, potential=espresso.interaction.LennardJones(epsilon, sigma, rc, shift))
  system.addInteraction(interaction)

  integrator     = espresso.integrator.VelocityVerlet(system)  
  integrator.dt  = dt
  if (temperature != None):
    thermostat             = espresso.integrator.LangevinThermostat(system)
    thermostat.gamma       = 1.0
    thermostat.temperature = temperature
    integrator.addExtension(thermostat)
  mass = 1.0
  if xyzrfilename: 
    new_particles = []
    props     = ['id', 'type', 'mass', 'pos', 'v', 'radius']
    for idx in range(num_particles):
      part = [ pidf[idx], typef[idx], mass,
               espresso.Real3D(xposf[idx],yposf[idx],zposf[idx]),
               espresso.Real3D(xvelf[idx],yvelf[idx],zvelf[idx]),
               radiusf[idx] ]
      new_particles.append(part)
      if idx % 1000 == 0:
        system.storage.addParticles(new_particles, *props)
        system.storage.decompose()
        new_particles = []      
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
  elif xyzfilename: 
    new_particles = []
    props     = ['id', 'type', 'mass', 'pos', 'v']
    for idx in range(num_particles):
      part = [ pidf[idx], typef[idx], mass,
               espresso.Real3D(xposf[idx],yposf[idx],zposf[idx]),
               espresso.Real3D(xvelf[idx],yvelf[idx],zvelf[idx])]
      new_particles.append(part)
      if idx % 1000 == 0:
        system.storage.addParticles(new_particles, *props)
        system.storage.decompose()
        new_particles = []      
    system.storage.addParticles(new_particles, *props)
    system.storage.decompose()
  else:  
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
