from cpp.espresso.world import _Particle
from wrap_cpp import choose

class Particle(_Particle):
    'This describes a particle.'
    def set(self, id=None, x=None, y=None, z=None) :
        'Set the parameters of a Particle.'
        _Particle.set(self, 
                      choose(id, self.id), 
                      choose(x, self.x),
                      choose(y, self.y),
                      choose(z, self.z))
        
    def __init__(self, id=None, x=None, y=None, z=None) :
        'Initialize a Particle.'
        _Particle.__init__(self)
        self.set(id, x, y, z)

    # Testing boost::optional template
    def set_optional(self, id=None, x=None, y=None, z=None) :
        print "Calling set_optional..."
        _Particle.set_optional(self, id, x, y, z)
