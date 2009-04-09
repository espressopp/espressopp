#######################################################
#
#   Test program for testing Boost-Python interface to 
#   the C++ class espresso::particles::Storage|Computer

#   particles_Storage        : espresso::particles::storage
#   particles_ParticleWriter : new added class in bindings for test here
#                              is an extension of expresso::particles::computer
#

import espresso
from _espresso import particles_Storage
from _espresso import particles_PythonComputer
from _espresso import Real3D
from _espresso import Real3DProperty
import random

class ParticleWriter(particles_PythonComputer):
    def __init__(self, _property, storage):
        self.property = _property
        particles_PythonComputer.__init__(self, storage)

    def each(self, pid):
        print("%d %s" % (pid, self.property[pid]))

storage = particles_Storage()

pos = Real3DProperty(storage)

for count in range(0,100):
    pid = storage.addParticle()
    pos[pid] = Real3D(random.random(), random.random(), random.random())

storage.foreach(ParticleWriter(pos, storage))

