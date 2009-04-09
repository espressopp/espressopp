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
from _espresso import IntegerProperty
import random

# tag particles in a small sphere around (0,0,0)
class ParticleTester(particles_PythonComputer):
    def __init__(self, _position, _tag, storage):
        particles_PythonComputer.__init__(self, storage)
        self.position = _position
        self.tag = _tag

    def each(self, pid):
        pos = self.position[pid]
        if pos*pos < 0.5:
            self.tag[pid] = 1
        else:
            self.tag[pid] = 0

class ParticleWriter(particles_PythonComputer):
    def __init__(self, _property, _tag, storage):
        particles_PythonComputer.__init__(self, storage)
        self.property = _property
        self.tag = _tag

    def each(self, pid):
        if self.tag[pid]:
            print("%d %s" % (pid, self.property[pid]))

storage = particles_Storage()

pos = Real3DProperty(storage)

for count in range(0,100):
    pid = storage.addParticle()
    pos[pid] = Real3D(random.random(), random.random(), random.random())

# add property a posteriori
tag = IntegerProperty(storage)

#tag particles
storage.foreach(ParticleTester(pos, tag, storage))

# and print tagged ones
storage.foreach(ParticleWriter(pos, tag, storage))

