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

#read in the restart configuration
pid, type, x, y, z, vx, vy, vz, Lx, Ly, Lz = espressopp.tools.readxyz('dump.86500.xyz')
num_particles 		= len(pid)
box			= (Lx, Ly, Lz)
num_chains		= 328
monomers_per_chain	= 10
rc 			= 2 * pow(2, 1./6.)
skin			= 0.3
epsilon			= 1.
sigma			= 1.
temperature		= 1.0
dt			= 0.005

# create default Lennard Jones (WCA) system with 0 particles and cubic box (L=40)
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
#thermostat     = espressopp.integrator.LangevinThermostat(system)
#thermostat.gamma  = 1.0
#thermostat.temperature = temperature
#integrator.addExtension(thermostat)

print 'timestep is ', integrator.dt
#print 'gamma of the thermostat is ', thermostat.gamma
#print 'temperature of the thermostat is ', thermostat.temperature

# redefine bonds
props    = ['id', 'type', 'mass', 'pos', 'v']
bondlist = espressopp.FixedPairList(system.storage)
mass = 1.0

for i in range(num_chains):
	chain = []
	bonds = []
	for k in range(monomers_per_chain):
		partid = i * monomers_per_chain + k + 1
		pos = Real3D(x[partid-1], y[partid-1], z[partid-1])
		v = Real3D(vx[partid-1], vy[partid-1], vz[partid-1])
		particle = [partid, type[partid-1], mass, pos, v]
		chain.append(particle)
		if partid % monomers_per_chain != 0:
			bonds.append((partid,partid + 1))
	system.storage.addParticles(chain, *props)
	system.storage.decompose()
	bondlist.addBonds(bonds)

system.storage.decompose()

potFENE   = espressopp.interaction.FENE(K=30.0, r0=0.0, rMax=1.5)
interFENE = espressopp.interaction.FixedPairListFENE(system, bondlist, potFENE)
system.addInteraction(interFENE)

for k in range(5):
	espressopp.tools.analyse.info(system, integrator)
	cmvel = Real3D(0.,0.,0.)
	for i in range(1,num_particles+1):
		particle = system.storage.getParticle(i)
		cmvel += particle.v
	print "centre of mass velocity is ", cmvel
	integrator.run(100)

#thermostat.disconnect() # disconnect md-thermostat as we want to run lb-md coupled system

# define a LB grid
lb = espressopp.integrator.LatticeBoltzmann(system, Ni=Int3D(16, 16, 16))
initPop = espressopp.integrator.LBInitPopUniform(system,lb)
initPop.createDenVel(1.0, Real3D(0.,0.,0.0))

# declare gammas responsible for viscosities (if they differ from 0)
lb.gamma_b = 0.5
lb.gamma_s = 0.5

# specify desired temperature (set the fluctuations if any)
lb.lbTemp = 0.000025
lb.fricCoeff = 20.

# add extension to the integrator
integrator.addExtension(lb)

print integrator.dt
print integrator.step
print thermostat.gamma
print thermostat.temperature
print lb.fricCoeff

#plt.figure()
T   = espressopp.analysis.Temperature(system)
#x   = []
#yT  = []
#yTmin = 0.2
#yTmax = 1.8

#plt.subplot(211)
#gT, = plt.plot(x, yT, 'ro')

lb.nSteps=5

# write output to a datafile
f = open('temp_L16_N328_G20_Nmd5_dt0.005.dat', 'a')

for k in range(500):
	lb.readCouplForces()
	integrator.run(100)
	currT = T.compute()
	s = str(integrator.step)
	f.write(s+'\t')
	mdoutput = 'dump.' + s + '.xyz'
	s = str(currT)
	f.write(s+'\n')
	lb.saveCouplForces()
	espressopp.tools.writexyz(mdoutput, system)
#	x.append(integrator.dt * integrator.step)
#	currT = T.compute()
#	yT.append(currT)
#	s = str(integrator.step)
#	f.write(s+'\t')
#	s = str(currT)
#	f.write(s+'\n')
#	plt.subplot(211)
#	plt.axis([x[0], x[-1], yTmin, yTmax ])
#	gT.set_ydata(yT)
#	gT.set_xdata(x)
#	plt.draw()
#
#
#plt.savefig('lb1.0_c1.0_L16_N328_G20_2.pdf')
f.close()
