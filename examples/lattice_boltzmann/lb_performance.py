# DEMONSTRATION OF THE LATTICE-BOLTZMANN SIMULATION
#
s=20
import espressopp
from espressopp import Int3D
from espressopp import Real3D
# create default Lennard Jones (WCA) system with 0 particles and cubic box (L=40)
system, integrator = espressopp.standard_system.LennardJones(100, box=(s, s, s), temperature=1.)
system.rng.seed = 123456
integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = 0.01

# define a LB grid
nodeGrid=espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
# add extension to the integrator
integrator.addExtension(lb)

lb.nSteps = 10
lb.visc_b = 3.
lb.visc_s = 3.
lb.profStep = 5000

#initPop = espressopp.integrator.LBInitPopUniform(system,lb)
initPop = espressopp.integrator.LBInitPopWave(system,lb)
initPop.createDenVel(1.0, Real3D(0.,0.,0.))
#initPop.createDenVel(1.0, Real3D(0.,0.,0.001))

lboutputScreen = espressopp.analysis.LBOutputScreen(system,lb)
OUT3=espressopp.integrator.ExtAnalyze(lboutputScreen,5000)
integrator.addExtension(OUT3)

# specify desired temperature (set the fluctuations if any)
lb.lbTemp = 1.
lb.nSteps = 1
## add some profiling statistics for the run
for k in range (3):
#	lb.readCouplForces()
	integrator.run(50000)
	s = str(integrator.step)
	mdoutput = 'dump.' + s + '.xyz'
	espressopp.tools.writexyz(mdoutput, system)
	lb.saveCouplForces()
