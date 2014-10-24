# DEMONSTRATION OF THE LATTICE-BOLTZMANN SIMULATION
#
import espresso
import cProfile, pstats
from espresso import Int3D
from espresso import Real3D

# create default Lennard Jones (WCA) system with 0 particles and cubic box (L=40)
system, integrator = espresso.standard_system.LennardJones(100, box=(20, 20, 20), temperature=1.2)

#integrator.dt = 0.000001
#integrator.run(10000)
#print "Finished with warming up"
#integrator.dt = 0.005

# define a LB grid
lb = espresso.integrator.LatticeBoltzmann(system, Ni=Int3D(20, 20, 20))
initPop = espresso.integrator.LBInitPopUniform(system,lb)
#initPop = espresso.integrator.LBInitPopWave(system,lb)
initPop.createDenVel(1.0, Real3D(0.,0.,0.0))

# declare gammas responsible for viscosities (if they differ from 0)
lb.gamma_b = 0.5
lb.gamma_s = 0.5

# specify desired temperature (set the fluctuations if any)
#lb.lbTemp = 0.0
lb.lbTemp = 0.000025

# add extension to the integrator
integrator.addExtension(lb)

# output velocity profile vz (x)
lboutputVzOfX = espresso.analysis.LBOutputProfileVzOfX(system,lb)
OUT1=espresso.integrator.ExtAnalyze(lboutputVzOfX,100)
integrator.addExtension(OUT1)

# output velocity vz at a certain lattice site as a function of time
#lboutputVzInTime = espresso.analysis.LBOutputVzInTime(system,lb)
#OUT2=espresso.integrator.ExtAnalyze(lboutputVzInTime,100)
#integrator.addExtension(OUT2)

# output onto the screen
lboutputScreen = espresso.analysis.LBOutputScreen(system,lb)
OUT3=espresso.integrator.ExtAnalyze(lboutputScreen,100)
integrator.addExtension(OUT3)

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
## add some profiling statistics for the run
cProfile.run("integrator.run(10000)",'profiler_stats')
p = pstats.Stats('profiler_stats')
p.strip_dirs().sort_stats("time").print_stats(10)
