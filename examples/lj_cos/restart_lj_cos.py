# DEMONSTRATION OF THE LJ-COS POTENTIAL CONTROLLING SOLVENT QUALITY
#
import espressopp
from espressopp import Int3D
from espressopp import Real3D

# initial parameters of the simulation read in from the "old" file
pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz('dump.150000.xyz')
box   			= (Lx, Ly, Lz)
num_chains 		= 70
mon_per_chain 	= 50
rc					= 1.5
skin				= 0.3
dt					= 0.005
sigma				= 1.
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

# redefine bonds
props    = ['id', 'type', 'mass', 'pos', 'v']
bondlist = espressopp.FixedPairList(system.storage)
mass     = 1.0

# SEPARATE BONDED AND NON-BONDED LJs
#exclusionlist = []
for i in range(num_chains):
	chain = []
	bonds = []
	for k in range(mon_per_chain):
		partid = i * mon_per_chain + k + 1
		pos = Real3D(x[partid-1], y[partid-1], z[partid-1])
		v = Real3D(vx[partid-1], vy[partid-1], vz[partid-1])
		particle = [partid, type[partid-1], mass, pos, v]
		chain.append(particle)
		if partid % mon_per_chain != 0:
			bonds.append((partid,partid + 1))
	system.storage.addParticles(chain, *props)
	system.storage.decompose()
	bondlist.addBonds(bonds)
# SEPARATE BONDED AND NON-BONDED LJs
#		if k<mon_per_chain-1: 
#			exclusionlist.append(bonds[k])
	
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

T   = espressopp.analysis.Temperature(system)

# write output to a datafile
f_temp = open('temperature_rest.dat', 'a')
f_stat = open('static_prop_rest.dat', 'a')

num_particles = num_chains * mon_per_chain

#integrator.step = 0

for k in range(150):
	integrator.run(1000)
	currT = T.compute()
	s = str(integrator.step)
	f_temp.write(s+'\t')
	f_stat.write(s+'\t')
	mdoutput = 'dump.' + s + '.xyz'
	espressopp.tools.writexyz(mdoutput, system, unfolded = True)
	s = str(currT)
	f_temp.write(s+'\n')
	espressopp.tools.analyse.info(system, integrator)

	# calculate end-to-end
	R2ee = 0.
	for i in range(num_chains):
		im_x = 0
		im_y = 0
		im_z = 0
		for j in range(mon_per_chain-1):
			p1 = system.storage.getParticle(i*mon_per_chain+j+1)
			p2 = system.storage.getParticle(i*mon_per_chain+j+2)
			dr = Real3D(p2.pos.x-p1.pos.x, p2.pos.y-p1.pos.y, p2.pos.z-p1.pos.z)
			if (dr.x < -.5 * Lx):
				im_x += 1
			elif (dr.x > .5 * Lx):
				im_x -= 1

			if (dr.y < -.5 * Ly):
				im_y += 1
			elif (dr.y > .5 * Ly):
				im_y -= 1

			if (dr.z < -.5 * Lz):
				im_z += 1
			elif (dr.z > .5 * Lz):
				im_z -= 1
			
		p1 = system.storage.getParticle(i*mon_per_chain+1)
		p2 = system.storage.getParticle((i+1)*mon_per_chain)
		Lbox = system.bc.boxL
		c1x = p1.pos.x
		c1y = p1.pos.y
		c1z = p1.pos.z
		c2x = p2.pos.x + Lx * im_x
		c2y = p2.pos.y + Ly * im_y
		c2z = p2.pos.z + Lz * im_z
		dr = Real3D(c2x-c1x, c2y-c1y, c2z-c1z)
		d2r = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]
		R2ee += d2r
	R2ee /= num_chains
	s = str(R2ee)
	f_stat.write(s+'\n')
	print R2ee
f_temp.close()
f_stat.close()
