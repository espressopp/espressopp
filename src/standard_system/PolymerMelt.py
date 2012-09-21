import espresso
import MPI

def PolymerMelt(num_chains, monomers_per_chain, box, bondlen=0.97, rc=1.12246, skin=0.3, dt=0.005, epsilon=1.0, sigma=1.0, shift='auto', temperature=None):
  '''
  returns random walk polymer melt system and integrator:
  if tempearture is != None then Langevin thermostat is set to temperature (gamma is 1.0)
  '''
  system         = espresso.System()
  system.rng     = espresso.esutil.RNG()
  system.bc      = espresso.bc.OrthorhombicBC(system.rng, box)
  system.skin    = skin
  nodeGrid       = espresso.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
  cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
  system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
  interaction    = espresso.interaction.VerletListLennardJones(espresso.VerletList(system, cutoff=rc+skin))
  interaction.setPotential(type1=0, type2=0, potential=espresso.interaction.LennardJones(epsilon, sigma, rc, shift))
  system.addInteraction(interaction)
  
  integrator     = espresso.integrator.VelocityVerlet(system)  
  integrator.dt  = dt
  if (temperature != None):
    thermostat             = espresso.integrator.LangevinThermostat(system)
    thermostat.gamma       = 1.0
    thermostat.temperature = temperature
    integrator.addExtension(thermostat)
  
  props    = ['id', 'type', 'mass', 'pos', 'v']
  vel_zero = espresso.Real3D(0.0, 0.0, 0.0)
  bondlist = espresso.FixedPairList(system.storage)
  pid      = 1
  type     = 0
  mass     = 1.0  
  chain    = []
  for i in range(num_chains):
    startpos = system.bc.getRandomPos()
    positions, bonds = espresso.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen)
    for k in range(monomers_per_chain):  
      part = [pid + k, type, mass, positions[k], vel_zero]
      chain.append(part)
    pid += monomers_per_chain
    type += 1
    system.storage.addParticles(chain, *props)
    system.storage.decompose()
    chain = []
    bondlist.addBonds(bonds)

  system.storage.decompose()

  # FENE bonds
  potFENE   = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
  interFENE = espresso.interaction.FixedPairListFENE(system, bondlist, potFENE)
  system.addInteraction(interFENE)
    
  return system, integrator
