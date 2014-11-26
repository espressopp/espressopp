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
#num_chains		= 328
num_chains		= 1
monomers_per_chain	= 10
L			= 16
box			= (L, L, L)
bondlen			= 0.97
rc 			= 2 * pow(2, 1./6.)
skin			= 0.3
dt			= 0.000001
epsilon			= 0.
sigma			= 1.
temperature		= 1.2
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

print integrator.dt
print thermostat.gamma
print thermostat.temperature

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

#force_capping   = espresso.integrator.CapForce(system, 500.0)
#integrator.addExtension(force_capping)
espresso.tools.analyse.info(system, integrator)

print "First phase of the warm up. Epsilon will be increased from 0. to 0.25"
new_epsilon = 0.
for k in range(10):
	integrator.run(2000)
	espresso.tools.analyse.info(system, integrator)
	new_epsilon += 0.025
	potLJ = espresso.interaction.LennardJones(new_epsilon, sigma, rc)
	interaction.setPotential(type1=0, type2=0, potential=potLJ)

print "Second phase of the warm up with timestep 0.0001. Epsilon will be increased from to 1."
print new_epsilon
integrator.dt = 0.00005
for k in range(10):
	integrator.run(1000)
	espresso.tools.analyse.info(system, integrator)
	new_epsilon += 0.075
	potLJ = espresso.interaction.LennardJones(new_epsilon, sigma, rc)
	interaction.setPotential(type1=0, type2=0, potential=potLJ)

print "Third phase of the warm up with timestep 0.005. Epsilon is 1."
print new_epsilon
integrator.dt = 0.005
for k in range(10):
	integrator.run(1000)
	espresso.tools.analyse.info(system, integrator)

thermostat.gamma  = 0.0
thermostat.temperature = 0.0

# define a LB grid
lb = espresso.integrator.LatticeBoltzmann(system, Ni=Int3D(16, 16, 16))
initPop = espresso.integrator.LBInitPopUniform(system,lb)
#initPop = espresso.integrator.LBInitPopWave(system,lb)
initPop.createDenVel(1.0, Real3D(0.,0.,0.0))

# declare gammas responsible for viscosities (if they differ from 0)
lb.gamma_b = 0.5
lb.gamma_s = 0.5

# specify desired temperature (set the fluctuations if any)
#lb.lbTemp = 0.0
lb.lbTemp = 0.000025
lb.fricCoeff = 20.
#lb.fricCoeff = 0.

# add extension to the integrator
integrator.addExtension(lb)

# output velocity profile vz (x)
#lboutputVzOfX = espresso.analysis.LBOutputProfileVzOfX(system,lb)
#OUT1=espresso.integrator.ExtAnalyze(lboutputVzOfX,100)
#integrator.addExtension(OUT1)

# output velocity vz at a certain lattice site as a function of time
#lboutputVzInTime = espresso.analysis.LBOutputVzInTime(system,lb)
#OUT2=espresso.integrator.ExtAnalyze(lboutputVzInTime,100)
#integrator.addExtension(OUT2)

# output onto the screen
#lboutputScreen = espresso.analysis.LBOutputScreen(system,lb)
#OUT3=espresso.integrator.ExtAnalyze(lboutputScreen,100)
#integrator.addExtension(OUT3)

print integrator.dt
print thermostat.gamma
print thermostat.temperature
print lb.fricCoeff
# set external constant (gravity-like) force
#lbforce = espresso.integrator.LBInitConstForce(system,lb)
#lbforce.setForce(Real3D(0.,0.,0.0001))
# run 500 steps with it
#integrator.run(500)
#integrator.run(100000)

# add a periodic force with a specified amplitude to the existing body force
#lbforce2 = espresso.integrator.LBInitPeriodicForce(system,lb)
#lbforce2.addForce(Real3D(0.,0.,0.0005))
#lb.lbTemp = 0.0000005
## run 500 steps with it
#integrator.run(500)
##
#plt.figure()
#T   = espresso.analysis.Temperature(system)
#x   = []
#yT  = []
#yTmin = 0.2
##x.append(integrator.dt * integrator.step)
##yT.append(T.compute())
#yTmax = 1.8

#plt.subplot(211)
#gT, = plt.plot(x, yT, 'ro')

# write output to a datafile
#f = open('temp_lb1.0_c1.0_L16_N328_G50.dat', 'a')

integrator.run(50)
#print "new"
integrator.run(50)
#for k in range(100):
#	print "before"
#	integrator.run(100)
#	print "after"
#	x.append(integrator.dt * integrator.step)
#	currT = T.compute()
#	yT.append(currT)
#	s = str(integrator.step)
#	f.write(s+'\t')
#	s = str(currT)
#	f.write(s+'\n')
##	yTmax = max(yT)
#	plt.subplot(211)
#	plt.axis([x[0], x[-1], yTmin, yTmax ])
#	gT.set_ydata(yT)
#	gT.set_xdata(x)
#	plt.draw()


#plt.savefig('lb1.0_c1.0_L16_N328_G50.pdf')
#f.close()

## add some profiling statistics for the run
##cProfile.run("integrator.run(10000)",'profiler_stats')
##p = pstats.Stats('profiler_stats')
##p.strip_dirs().sort_stats("time").print_stats(10)
