# DEMONSTRATION OF THE LJ-COS POTENTIAL CONTROLLING SOLVENT QUALITY
#
import espressopp
from espressopp import Int3D
from espressopp import Real3D

# initial parameters of the simulation
num_chains 		= 16
mon_per_chain 	= 200
L 					= 15.55644475
#num_chains 		= 70
#mon_per_chain 	= 50
#L 					= 16.02813677
box 				= (L, L, L)
bondlen 			= 0.97
rc					= 1.5
skin				= 0.3
dt					= 0.0000001
sigma				= 0.
temperature		= 1.

system         = espressopp.System()
system.rng     = espressopp.esutil.RNG()
system.bc      = espressopp.bc.OrthorhombicBC(system.rng, box)
system.skin    = skin
nodeGrid       = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
cellGrid       = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

# set-up of the integrator and its timestep
integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = dt

# set-up of the thermostat
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
# SEPARATE BONDED AND NON-BONDED LJs
#exclusionlist = []
for i in range(num_chains):
	startpos = system.bc.getRandomPos()
	positions, bonds = espressopp.tools.topology.polymerRW(pid, startpos, mon_per_chain, bondlen)
	for k in range(mon_per_chain):
		part = [pid + k, type, mass, positions[k], vel_zero]
		chain.append(part)
# SEPARATE BONDED AND NON-BONDED LJs
#		if k<mon_per_chain-1: 
#			exclusionlist.append(bonds[k])
	pid += mon_per_chain
	system.storage.addParticles(chain, *props)
	system.storage.decompose()
	chain = []
	bondlist.addBonds(bonds)
	
# SEPARATE BONDED AND NON-BONDED LJs
#vl.exclude(exclusion_list)

system.storage.decompose()

# initial set-up of the Lennard-Jones cosine
vl      = espressopp.VerletList(system, cutoff = rc)
potLJcos = espressopp.interaction.LJcos(phi=0.0)
potLJcos.sigma = sigma
interLJcos = espressopp.interaction.VerletListLJcos(vl)
interLJcos.setPotential(type1=0, type2=0, potential=potLJcos)
system.addInteraction(interLJcos)

# FENE-potential
potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, bondlist, potFENE)
system.addInteraction(interFENE)

# force-capping
new_force_cap = 5000.
force_capping = espressopp.integrator.CapForce(system, new_force_cap)
integrator.addExtension(force_capping)
espressopp.tools.analyse.info(system, integrator)

print "First phase of the warm up. Sigma will be increased from 0. to 1.0 and timestep to 0.001"
new_sigma = 0.
dt_step = .25 * integrator.dt 
for l in range(4):
	print "start increasing sigma from ", new_sigma
	for k in range(10):
		for j in range (200):
			new_sigma += 0.000125
			potLJcos.sigma = new_sigma
			interLJcos.setPotential(type1=0, type2=0, potential=potLJcos)
			integrator.run(100)
		espressopp.tools.analyse.info(system, integrator)
	print "finished increasing sigma upto ", new_sigma
	print "start force de-capping from ", new_force_cap
	for k in range(40):
		integrator.dt += dt_step
		integrator.run(500)
		espressopp.tools.analyse.info(system, integrator)
	for k in range(40):
		new_force_cap += 100
		force_capping.setAbsCapForce(new_force_cap)
		integrator.run(500)
	print "finished force de-capping upto ", new_force_cap
	print "integrator timestep is ", integrator.dt

force_capping.disconnect()
print "switching off force capping"

print "new_sigma is ", new_sigma

print "Second phase of the warm up with a production timestep. Force capping is turned off."
integrator.dt = 0.005
for k in range(100):
	integrator.run(10000)
	espressopp.tools.analyse.info(system, integrator)

T   = espressopp.analysis.Temperature(system)

# write output to a datafile
f_temp = open('temperature.dat', 'a')
f_stat = open('static_prop.dat', 'a')

num_particles = num_chains * mon_per_chain

integrator.step = 0

for k in range(150):
	integrator.run(1000)
	currT = T.compute()
	s = str(integrator.step)
	f_temp.write(s+'\t')
	f_stat.write(s+'\t')
	mdoutput = 'dump.' + s + '.xyz'
	espressopp.tools.writexyz(mdoutput, system)
	s = str(currT)
	f_temp.write(s+'\n')
	espressopp.tools.analyse.info(system, integrator)

	# calculate end-to-end
	R2ee = 0.
	for i in range(num_chains):
		p1 = system.storage.getParticle(i*mon_per_chain+1)
		p2 = system.storage.getParticle((i+1)*mon_per_chain)
		Lbox = system.bc.boxL
		c1x = p1.pos.x + Lbox.x * p1.imageBox.x
		c1y = p1.pos.y + Lbox.y * p1.imageBox.y
		c1z = p1.pos.z + Lbox.z * p1.imageBox.z
		c2x = p2.pos.x + Lbox.x * p2.imageBox.x
		c2y = p2.pos.y + Lbox.y * p2.imageBox.y
		c2z = p2.pos.z + Lbox.z * p2.imageBox.z
		dr = Real3D(c1x-c2x, c1y-c2y, c1z-c2z)
		d2r = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]
		R2ee += d2r
	R2ee /= num_chains
	s = str(R2ee)
	f_stat.write(s+'\n')
	print R2ee
f_temp.close()
f_stat.close()
