### Test for Lennard Jones interaction
#
# ToDo JH: make a Python unit test of it

from interaction import *

from _espresso import Real3D

# use constructor with default values

lj = LennardJones()

# print the (default) values

print 'sigma = ', lj.getSigma()
print 'epsilon = ', lj.getEpsilon()
print 'cutoff = ', lj.getCutoff()

# compute some energies for distances in a loop

for r in (0.8, 0.9, 1.0, 1.1, 1.2, 1.3):
   print 'r = ', r, ', energy = ', lj.computeEnergy(r)

# Note: Real3D.tuple() returns a Python tuple with three float 

x = Real3D(1.3, 1.1, -0.8)

f = lj.computeForce(x)

print 'computeForce: x = ', x.tuple(), ', force = ', f.tuple()

