# Test for Lennard Jones interaction

from espresso.interaction import LennardJones
from espresso.base import Real3D
from espresso.integrator import VelocityVerlet

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

vvIntegrator = VelocityVerlet(0.01)
#vvIntegrator.setTimeStep(0.05)
