def choose(val, altval) :
    if (val is None) :
        return altval
    else :
        return val


# wrap LennardJones
from _espresso import interaction_LennardJones as _LennardJones 
class LennardJones(_LennardJones) :
    'The Lennard-Jones interaction.'
    super = _LennardJones

    # define setter
    def set(self, epsilon=None, sigma=None, cutoff=None) :
        # call the C++ setter
        self.super.set(self,
            choose(epsilon, self.epsilon),
            choose(sigma, self.sigma),
            choose(cutoff, self.cutoff)
            )

    # define single property setters
    # avoid using these if possible
    def setEpsilon(self, _epsilon) :
        set(epsilon=_epsilon)
    def setSigma(self, _sigma) : 
        set(sigma=_sigma)
    def setCutoff(self, _cutoff) : 
        set(cutoff=_cutoff)

    # define properties
    epsilon = property(super.getEpsilon, setEpsilon)
    sigma = property(super.getSigma, setSigma)
    cutoff = property(super.getCutoff, setCutoff)

    # define constuctor
    def __init__(self, epsilon=None, sigma=None, cutoff=None) :
        self.super.__init__(self)
        self.set(epsilon, sigma, cutoff)

    def computeEnergy(self, r) :
        'Compute and return the energy at the radius r.'
        return self.super.computeEnergy(self, r)

# wrap FENE
from _espresso import interaction_FENE as _FENE
class FENE(_FENE) :
    'The FENE interaction.'
    super = _FENE

    # define setter
    def set(self, K=None, r0=None, rMax=None) :
        # call the C++ setter
        self.super.set(self,
            choose(K, self.K),
            choose(r0, self.r0),
            choose(rMax, self.rMax)
            )

    # define single property setters
    # avoid using these if possible
    def setK(self, _K) :
        set(K=_K)
    def setR0(self, _r0) : 
        set(r0=_r0)
    def setRMax(self, _rMax) : 
        set(rMax=_rMax)

    # define properties
    K = property(super.getK, setK)
    r0 = property(super.getR0, setR0)
    rMax = property(super.getRMax, setRMax)

    # define constuctor
    def __init__(self, K=None, r0=None, rMax=None) :
        self.super.__init__(self)
        self.set(K, r0, rMax)

    def computeEnergy(self, r) :
        'Compute and return the energy at the radius r.'
        return self.super.computeEnergy(self, r)
    
