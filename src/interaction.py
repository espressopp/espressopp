# wrap LennardJones
from _espresso import interaction_LennardJones as _LennardJones 
class LennardJones(_LennardJones) :
    'The Lennard-Jones interaction.'
    super = _LennardJones

    cutoff = property(super.getCutoff, super.setCutoff)
    sigma = property(super.getSigma, super.setSigma)
    epsilon = property(super.getEpsilon, super.setEpsilon)

    def __init__(self, cutoff=None, sigma=None, epsilon=None) :
        self.super.__init__(self)
        self.set(cutoff=cutoff, sigma=sigma, epsilon=epsilon)

    def set(self, cutoff=None, sigma=None, epsilon=None) :
        if (cutoff is not None): self.cutoff = cutoff
        if (epsilon is not None): self.epsilon = epsilon
        if (sigma is not None): self.sigma = sigma

    def computeEnergy(self, r) :
        'Compute and return the energy at the radius r.'
        return self.super.computeEnergy(self, r)

# wrap FENE
from _espresso import interaction_FENE as _FENE
class FENE(_FENE) :
    'The FENE interaction.'
    super = _FENE

    K = property(super.getK, super.setK)
    r0 = property(super.getR0, super.setR0)
    rMax = property(super.getRMax, super.setRMax)

    def __init__(self, K=None, r0=None, rMax=None) :
        self.super.__init__(self)
        self.set(K=K, r0=r0, rMax=rMax)

    def set(self, K=None, r0=None, rMax=None) :
        if (K is not None): self.K = K
        if (r0 is not None): self.r0 = r0
        if (rMax is not None): self.rMax = rMax

    def computeEnergy(self, r) :
        'Compute and return the energy at the radius r.'
        return self.super.computeEnergy(self, r)
    
