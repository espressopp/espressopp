from espresso._espresso import interaction_LennardJones
_LennardJones=interaction_LennardJones

class LennardJones(_LennardJones):
    'The Lennard-Jones interaction.'
    def getCutoff(self) :
        'Returns the cutoff.'
        return super(LennardJones, self).getCutoff(self)
    def getSigma(self) :
        'Returns sigma.'
        return super(LennardJones, self).getSigma(self)
    def getEpsilon(self) :
        'Returns epsilon.'
        return super(LennardJones, self).getEpsilon(self)
