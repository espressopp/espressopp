# DEMONSTRATION OF THE LATTICE-BOLTZMANN SIMULATION
#
import espressopp
from espressopp import Int3D
from espressopp import Real3D

L = 20
tempT = 1.
num_particles = 100
# create default Lennard Jones (WCA) system with 0 particles and cubic box
system, integrator = espressopp.standard_system.LennardJones(num_particles, box=(L, L, L), temperature=tempT)
integrator     = espressopp.integrator.VelocityVerlet(system)
integrator.dt  = 0.005

# define a LB grid
nodeGrid=espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
integrator.addExtension(lb)

# set parameters of LB fluid (LJ units)
lb.lbTemp = tempT       # desired temperature
lb.nSteps = 10          # time step contrast between LB and MD (t_{LB} / t_{MD})
lb.visc_b = 3.          # bulk viscosity of LB fluid
lb.visc_s = 3.          # shear viscosity of LB fluid
lb.profStep = 5000      # time profiling frequency

# initialize populations
initDen = 1.
initVel = Real3D (0.)
#initPop = espressopp.integrator.LBInitPopWave(system,lb)   # sin-like in z-dir
initPop = espressopp.integrator.LBInitPopUniform(system,lb) # uniform
initPop.createDenVel( initDen, initVel )

# screen output
lboutputScreen = espressopp.analysis.LBOutputScreen(system,lb)
outScreen=espressopp.integrator.ExtAnalyze(lboutputScreen,lb.profStep)
integrator.addExtension(outScreen)

# run simulations
for k in range (3):
    integrator.run(50000)
    s = str(integrator.step)
    # output md configuration
    mdoutput = 'dump.' + s + '.xyz'
    espressopp.tools.fastwritexyz(mdoutput, system)
    # output LB configuration
    lb.keepLBDump()     # flag to keep previously saved LB state
    lb.saveLBConf()     # saves current state of the LB fluid
