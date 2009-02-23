from _espresso import interaction_LennardJones as _LennardJones

class LennardJones:
    'The Lennard-Jones interaction.'
    global _lj
    def __init__(self) :
        self._lj = _LennardJones()

    def getCutoff(self) :
        'Returns the cutoff.'
        return self._lj.getCutoff()

    def getSigma(self) :
        'Returns sigma.'
        return self._lj.getSigma()

    def getEpsilon(self) :
        'Returns epsilon.'
        return self._lj.getEpsilon()

    def computeEnergy(self, r) :
        'Compute and return the energy at the radius r.'
        return self._lj.computeEnergy(r)
