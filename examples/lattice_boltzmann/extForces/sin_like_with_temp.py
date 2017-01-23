# STOCHASTIC LATTICE-BOLTZMANN SIMULATION 
# WITH EXTERNAL GRAVITY- AND SIN-LIKE FORCES 
# AND FILE-OUTPUT OF THE VELOCITY COMPONENT VZ DEPENDENCE ON THE X COORDINATE

import espressopp
from espressopp import Int3D
from espressopp import Real3D

L = 20                          # box dimension in every direction
box = (L, L, L)
# create a default system with 0 particles and cubic box
system, integrator = espressopp.standard_system.Default(box=box)

# calculate CPU nodeGrid based on the number of CPUs
nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)

# LATTICE BOLTZMANN (LB) INITIALIZATION
# define an lb object and connect to the integrator
lb = espressopp.integrator.LatticeBoltzmann(system, nodeGrid)
integrator.addExtension(lb)

# set up temperature for thermal fluctuations ( LJ units )
lb.lbTemp = 1.

# initialize initial density and velocity
initDen = 1.                    # initial density on a site (LB units!)
initVel = Real3D(0.)            # initial velocity on a site (LB units!)
initPop = espressopp.integrator.LBInitPopUniform(system,lb)
initPop.createDenVel(initDen, initVel)

# output of the progress to the screen
outStep = 250                   # period of the statistics update
lbOutScreen    = espressopp.analysis.LBOutputScreen(system,lb)
outputToScreen = espressopp.integrator.ExtAnalyze(lbOutScreen, outStep)
integrator.addExtension(outputToScreen)

# output velocity profile vz (x)
outStep = 500
lbOutVelFile   = espressopp.analysis.LBOutputVzOfX(system,lb)
outputToFile   = espressopp.integrator.ExtAnalyze(lbOutVelFile, outStep)
integrator.addExtension(outputToFile)

# NOW WE WOULD LIKE TO APPLY TO LB LIQUID SOME EXTERNAL FORCES 
# set external constant (gravity-like) force
extForceToSet = Real3D(0., 0., 0.0001)
lbforce = espressopp.integrator.LBInitConstForce(system,lb)
lbforce.setForce( extForceToSet )

# run integrator for 500 steps with it
steps = 500
integrator.run( steps )

# NOW WE WOULD LIKE TO ADD SOME SIN-LIKE EXTERNAL FORCE TO THE EXISTING ONE 
# SIN MODULATION IS DONE ALONG Z-DIRECTION AS A FUNCTION OF X-COORDINATE
ampVz = 0.0001              # set amplitude of the sin like force in z direction
bodyFx = bodyFy = 0.        # body forces in x- and y-dirs to add (0. for now)
lbforceSin = espressopp.integrator.LBInitPeriodicForce(system, lb)
lbforceSin.addForce( Real3D(bodyFx, bodyFy, ampVz) )
integrator.run( steps )

# add an external sin-like force to the existing one
extForceToAdd = Real3D(0., 0., 0.0005)
lbforceSin.addForce( extForceToAdd )
integrator.run( steps )

# IF YOU LEAVE THE SYSTEM LIKE THIS THE VELOCITY OF THE FLUID WILL INCREASE
# WITH TIME AND EVERYTHING WILL EXPLODE EVENTUALLY. SO LETS PREVENT THIS
# DISASTER BY SETTING EXTERNAL FORCES TO ZERO AND RUNNING THE SYSTEM ONCE MORE
lbforce.setForce( Real3D( 0. ) )
integrator.run( steps + 1) # + 1 to print the last step profile

# THE WHOLE PROCESS CAN BE NOW SEEN IN THE FILES CONTAININT THE PROFILE VZ OF X
