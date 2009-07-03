from espresso.potential import CentralPotentialLocalBase

class MyPotential(CentralPotentialLocalBase):
    def __init__(self):
        CentralPotentialLocalBase.__init__(self)
    def computeEnergySqr(self, distSqr):
        return distSqr

mi=MyPotential()
print(mi.__class__)
mi.computeEnergy(2.0)
