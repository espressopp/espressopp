# DEMONSTRATION OF THE LATTICE-BOLTZMANN SIMULATION
#
import espresso
import cProfile, pstats
from espresso import Int3D
from espresso import Real3D
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

system         = espresso.System()
system.rng     = espresso.esutil.RNG()
system.bc      = espresso.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
nodeGrid       = espresso.tools.decomp.nodeGrid(espresso.MPI.COMM_WORLD.size)
cellGrid       = espresso.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espresso.storage.DomainDecomposition(system, nodeGrid, cellGrid)
interaction    = espresso.interaction.VerletListLennardJones(espresso.VerletList(system, cutoff=rc))
potLJ          = espresso.interaction.LennardJones(epsilon, sigma, rc)
interaction.setPotential(type1=0, type2=0, potential=potLJ)
system.addInteraction(interaction)

integrator     = espresso.integrator.VelocityVerlet(system)
integrator.dt  = dt
thermostat     = espresso.integrator.LangevinThermostat(system)
thermostat.gamma  = 1.0
thermostat.temperature = temperature
integrator.addExtension(thermostat)

print 'timestep is ', integrator.dt
print 'gamma of the thermostat is ', thermostat.gamma
print 'temperature of the thermostat is ', thermostat.temperature

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

potFENE   = espresso.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espresso.interaction.FixedPairListFENE(system, bondlist, potFENE)
system.addInteraction(interFENE)

force_capping = espresso.integrator.CapForce(system, 1000.0)
integrator.addExtension(force_capping)
espresso.tools.analyse.info(system, integrator)

print "First phase of the warm up. Epsilon will be increased from 0. to 1.0 and timestep to 0.001"
new_epsilon = 0.
new_force_cap = 10000.
for l in range(4):
	new_dt = integrator.dt
	for k in range(10):
		for j in range (100):
			integrator.run(10)
			new_epsilon += 0.00025
			potLJ = espresso.interaction.LennardJones(new_epsilon, sigma, rc)
			interaction.setPotential(type1=0, type2=0, potential=potLJ)
		espresso.tools.analyse.info(system, integrator)
		print "new_epsilon is ", new_epsilon, "sigma is ", sigma
	for k in range(9):
		integrator.run(1000)
		integrator.dt += new_dt
		espresso.tools.analyse.info(system, integrator)
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
	espresso.tools.analyse.info(system, integrator)

#thermostat.disconnect() # disconnect md-thermostat as we want to run lb-md coupled system

# define a LB grid
#lb = espresso.integrator.LatticeBoltzmann(system, Ni=Int3D(16, 16, 16))
#initPop = espresso.integrator.LBInitPopUniform(system,lb)
#initPop.createDenVel(1.0, Real3D(0.,0.,0.0))

# add extension to the integrator
#integrator.addExtension(lb)

T   = espresso.analysis.Temperature(system)

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
	espresso.tools.writexyz(mdoutput, system)
	s = str(currT)
	f.write(s+'\n')
	espresso.tools.analyse.info(system, integrator)
	cmvel = Real3D(0.,0.,0.)
	for i in range(1,num_particles+1):
		particle = system.storage.getParticle(i)
		cmvel += particle.v
	print "centre of mass velocity is ", cmvel
f.close()
