import pycpp

p=pycpp.Particle();

p.view()
p.set('hello',1,2,3)
p.view()
p.set(id='This is a test', posx=2.5);
p.view()
p.setid('nothing')
p.x=0.3
p.y=0.4
p.z=0.5
print p.x, p.y, p.z
