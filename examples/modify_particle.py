import espresso
from espresso import Real3D
from math import sin, cos

system, integrator = espresso.standard_system.Default((10, 10, 10))
dd = system.storage
dd.addParticle(0, Real3D(1.0, 1.0, 1.0))

dd.addParticles([ [1, Real3D(5.0, 5.0, 5.0)], \
                  [2, Real3D(5.0, 6.0, 7.0)], \
                  [3, Real3D(6.0, 1.0, 3.0)] ], 'id', 'pos')

dd.decompose()

sock=espresso.tools.vmd.connect(system)
print "start position:", dd.getParticle(0).pos
t = 0
x = dd.getParticle(0).pos[0]
y = dd.getParticle(0).pos[1]
z = dd.getParticle(0).pos[2]
while t < 10000:
  dd.modifyParticle(0,'pos', Real3D(x + t * 0.01, y + sin(0.02*t), z + cos(0.02*t)))
  espresso.tools.vmd.imd_positions(system, sock)
  t += 1
print "  end position:", dd.getParticle(0).pos
