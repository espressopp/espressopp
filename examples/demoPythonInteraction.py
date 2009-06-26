from espresso.interaction import CentralInteractionLocalBase

class MyInteraction(CentralInteractionLocalBase):
    def __init__(self):
        pass
    def computeEnergySqr(self, distSqr):
        return distSqr

mi=MyInteraction()
print(mi.__class__)
mi.computeEnergy(2.0)
