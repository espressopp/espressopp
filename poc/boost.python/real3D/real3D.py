import poc
import _real3D

class Real3D(_real3D.Real3D) :
    # string conversion
    def __str__(self) :
        return str(tuple(self))
    
    def __repr__(self) :
        return 'Real3D' + str(tuple(self))

r = Real3D(1.0, 2.0, 3.0)
print(r)

r[1] = 3.0
print(r)

_real3D.fold(r)
print(r)
