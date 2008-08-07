# LD_LIBRARY_PATH=$(BOOST_LIB)
import sys
sys.path.append('python')

import espresso.world

p = espresso.world.Particle(y=2.0)
print "Particle: ", p   

p.set(x=3.0,id='BLA')
print "Particle: ", p   

#p.set_optional('BLUB', 1.0, 2.0, 3.0)
#print "Particle: ", p

#p.set_optional(id='BLABER', x=1.0, z=3.0)
#print "Particle: ", p   

