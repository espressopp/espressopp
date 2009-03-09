#######################################################
#
#   Test program for testing Boost-Python interface to 
#   the C++ class espresso::particles::Storage|Computer

#   particles_Storage        : espresso::particles::storage
#   particles_ParticleWriter : new added class in bindings for test here
#                              is an extension of expresso::particles::computer
#

from _espresso import particles_Storage
from _espresso import particles_ParticleWriter

# ToDo: is it possible to make something like this
# class MyWriter(particles_Computer)
#    and overwrite () operator

storage = particles_Storage()

# Be careful: the default argument size_t dim = 1 does not export to Python

pos = storage.addPropertyReal3D(1)

# we use a routine of Storage to fill up with particles

size = 3.3
N    = 5

storage.fillWithLattice(size, N, pos)

# make a object of the test class that will print positions

writer = particles_ParticleWriter(storage, pos)

storage.foreach(writer)

