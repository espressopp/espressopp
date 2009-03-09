from _espresso import esutil_Real3D 

class Real3D(esutil_Real3D) :
    def tuple(self) :
        return (self[0], self[1], self[2])

    def __str__(self) :
        return str(self.tuple())
