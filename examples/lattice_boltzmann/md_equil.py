# DEMONSTRATION OF THE LATTICE-BOLTZMANN SIMULATION
#
import espressopp
import cProfile, pstats
from espressopp import Int3D
from espressopp import Real3D
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

# create default Lennard Jones (WCA) system with 0 particles and cubic box (L=40)
num_chains		= 328
monomers_per_chain	= 10
L			= 16
box			= (L, L, L)
bondlen			= 0.97
rc 			= 2 * pow(2, 1./6.)
skin			= 0.3
dt			= 0.0000001
epsilon			= 0.
sigma			= 1.
temperature		= 1.0
print "Initial values"

system         = espressopp.System()
system.rng     = espressopp.esutil.RNG()
system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
nodeGrid       = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)
interaction    = espressopp.interaction.VerletListLennardJones(espressopp.VerletList(system, cutoff=rc))
potLJ          = espressopp.interaction.LennardJones(epsilon, sigma, rc)
interaction.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interaction)

integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = dt
thermostat     = espressopp.integrator.LangevinThermostat(system)
thermostat.gamma  = 1.0
thermostat.temperature = temperature
integrator.addExtension(thermostat)

print 'timestep is ', integrator.dt
print 'gamma of the thermostat is ', thermostat.gamma
print 'temperature of the thermostat is ', thermostat.temperature

props    = ['id', 'type', 'mass', 'pos', 'v']
vel_zero = espressopp.Real3D(0.0, 0.0, 0.0)

bondlist = espressopp.FixedPairList(system.storage)
pid      = 1
type     = 0
mass     = 1.0
chain    = []

for i in range(num_chains):
	startpos = system.bc.getRandomPos()
	positions, bonds = espressopp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen)
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

potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, bondlist, potFENE)
system.addInteraction(interFENE)

force_capping = espressopp.integrator.CapForce(system, 1000.0)
integrator.addExtension(force_capping)
espressopp.tools.analyse.info(system, integrator)

print "First phase of the warm up. Epsilon will be increased from 0. to 1.0 and timestep to 0.001"
new_epsilon = 0.
new_force_cap = 10000.
for l in range(4):
	new_dt = integrator.dt
	for k in range(10):
		for j in range (100):
			integrator.run(10)
			new_epsilon += 0.00025
			potLJ = espressopp.interaction.LennardJones(new_epsilon, sigma, rc)
			interaction.setPotential(type1=0, type2=0, potential=potLJ)
		espressopp.tools.analyse.info(system, integrator)
		print "new_epsilon is ", new_epsilon, "sigma is ", sigma
	for k in range(9):
		integrator.run(1000)
		integrator.dt += new_dt
		espressopp.tools.analyse.info(system, integrator)
		new_force_cap += 1000
		force_capping.setAbsCapForce(new_force_cap)
		print "new force capping is ", new_force_cap
	print "integrator timestep is ", integrator.dt

force_capping.disconnect()
print "switching off force capping"

print "new_epsilon is ", new_epsilon, "sigma is ", sigma

print "Second phase of the warm up with a production timestep. Force capping is turned off."
integrator.dt = 0.005
for k in range(10):
	integrator.run(1000)
	espressopp.tools.analyse.info(system, integrator)

#thermostat.disconnect() # disconnect md-thermostat as we want to run lb-md coupled system

# define a LB grid
#lb = espressopp.integrator.LatticeBoltzmann(system, Ni=Int3D(16, 16, 16))
#initPop = espressopp.integrator.LBInitPopUniform(system,lb)
#initPop.createDenVel(1.0, Real3D(0.,0.,0.0))

# add extension to the integrator
#integrator.addExtension(lb)

T   = espressopp.analysis.Temperature(system)

# write output to a datafile
f = open('temperature.dat', 'a')

num_particles = num_chains * monomers_per_chain

#for k in range(500):
for k in range(5):
	integrator.run(100)
	currT = T.compute()
	s = str(integrator.step)
	f.write(s+'\t')
	mdoutput = 'dump.' + s + '.xyz'
	espressopp.tools.writexyz(mdoutput, system)
	s = str(currT)
	f.write(s+'\n')
	espressopp.tools.analyse.info(system, integrator)
	cmvel = Real3D(0.,0.,0.)
	for i in range(1,num_particles+1):
		particle = system.storage.getParticle(i)
		cmvel += particle.v
	print "centre of mass velocity is ", cmvel
f.close()
