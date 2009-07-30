from espresso.potential import CentralPotentialLocal

if __name__ == 'espresso.pmi':
    class MyPotentialLocal(CentralPotentialLocal):
        def __init__(self, a, b):
            CentralPotentialLocalBase.__init__(self)
            self.a = a
            self.b = b

        def computeEnergySqr(self, distSqr):
            return self.a*distSqr - self.b*distSqr*distSqr

        def computeForce(self, dist):
            return self.a * dist

else:

    .
    .
    allpairs = ...
    myint = potential.Interaction(
        potential=pmi.create('MyPotentialLocal', a, b),
        pairs=allpairs)
    myint.connect(integrator)
                          

