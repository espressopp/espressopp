#  Python example script that does the same as forceloop but with
#  epxorted routines from the C++ classes
#
#  Authors: Jon + Thomas
#
#  ToDo's: 
#     - bindings for particle/All and particel/Set (JH)
#     - addtional bindings for integrator/VerlocityVerlet (JH)
#     - ParticleWriter as PythonComputer (TB)
#     - bindings for forcecompute (TB)
#
#  Furter ToDo's (after this script running)
#     - Extend all C++ classes in xxx/__init__.py as Local + PMI classes
#     - writing UnitTest in Python

import random
import espresso

from _espresso import particles_Storage as Storage
from _espresso import Real3DProperty
from _espresso import Real3D
from _espresso import bc_PBC as PBC
# from _espresso import particles_All as All
from _espresso import pairs_All as All2
from _espresso import interaction_LennardJones as LennardJones


particleStorage = Storage()
position = Real3DProperty(particleStorage)
velocity = Real3DProperty(particleStorage)
force    = Real3DProperty(particleStorage)

N = 3
SIZE = 4.0

for i in range(N):
   for j in range(N):
      for k in range(N):
     
       r = 0.4 + 02 * random.random()
       x = (i + r) / N * SIZE;
       y = (j + r) / N * SIZE; 
       z = (k + r) / N * SIZE;

       id = particleStorage.addParticle()

       position[id] = Real3D(x, y, z);
       velocity[id] = Real3D(x, y, z);
       force[id] = Real3D(0.0)

pbc = PBC()
pbc.set(SIZE)

allSet = All(Storage)

allpairs = All2(pbc, allSet, position)

ljint = LennardJones();

ljint.set(1.0, 1.0, 2.5)

# not yet
# forcecompute = ljint->createForceComputer(ForceComputer(forceRef))

allpairs.foreach(forcecompute)

particleStorage.foreach(pWriter);

integrator = VelocityVerlet(allSet, position, velocity, force)

integrator.setTimeStep(0.005)
integrator.addForce(ljint, allpairs)
integrator.run(100)

particleStorage.foreach(pWriter)

