### Test for Lennard Jones interaction
#
# ToDo JH: make a Python unit test of it

from espresso.interaction import LennardJones
from espresso.esutil import Real3D

# use constructor with default values

lj = LennardJones()

# print the (default) values

print('sigma = %f' % lj.sigma)
print('epsilon = %f' % lj.epsilon)
print('cutoff = %f' % lj.cutoff)

# compute some energies for distances in a loop
for r in (0.8, 0.9, 1.0, 1.1, 1.2, 1.3):
   print('r = %f, energy = %f' % (r, lj.computeEnergy(r)))

x = Real3D(1.3, 1.1, -0.8)
f = lj.computeForce(x)

print('computeForce: x = %s, force = %s' % (x, f))

