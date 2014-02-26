# DEMONSTRATION OF THE LATTICE-BOLTZMANN SIMULATION
import espresso
from espresso import Int3D
from espresso import Real3D

# create a default system with 0 particles and cubic box
system, integrator = espresso.standard_system.Default(box=(20, 20, 20))

# LATTICE BOLTZMANN INITIALIZATION
# define grid and connect to the integrator
lb = espresso.integrator.LatticeBoltzmann(system, Ni=Int3D(20, 20, 20))
integrator.addExtension(lb)

# output of the progress to the screen
lboutputScreen = espresso.analysis.LBOutputScreen(system,lb)
OUT1=espresso.integrator.ExtAnalyze(lboutputScreen,100)
integrator.addExtension(OUT1)

# initialize initial density and velocity
initPop = espresso.integrator.LBInitPopUniform(system,lb)
initPop.createDenVel(1.0, Real3D(0.,0.,0.))

# APPLICATION OF FORCES TO THE LIQUID
# output velocity profile vz (x)
lboutputVzOfX = espresso.analysis.LBOutputProfileVzOfX(system,lb)
OUT2=espresso.integrator.ExtAnalyze(lboutputVzOfX,100)
integrator.addExtension(OUT2)

# set external constant (gravity-like) force
lbforce = espresso.integrator.LBInitConstForce(system,lb)
lbforce.setForce(Real3D(0.,0.,0.0001))
# run 500 steps with it
integrator.run(500)

# set external constant (gravity-like) force to zero
lbforce.addForce(Real3D(0.,0.,0.00005))
# run 500 steps with it
integrator.run(500)
