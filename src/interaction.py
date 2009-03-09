import pmi
import espresso.esutil

# pmi import _espresso
pmi.exec_('from _espresso import interaction_LennardJones')

def choose(val, altval) :
    if (val is None) :
        return altval
    else :
        return val


# wrap LennardJones
class LennardJones (object):
    'The Lennard-Jones interaction.'

    def __init__(self, epsilon = 1.0, sigma = 1.0, cutoff = 2.0) :
        """Initialize the (parallel) Lennard Jones object. 

        The parameters are identical to set."""
        object.__init__(self)
        # create the pmi object
        self.worker = pmi.create('interaction_LennardJones')
        # set the defaults
        self.set(epsilon, sigma, cutoff)

    # define setter
    def set(self, epsilon=None, sigma=None, cutoff=None) :
        'Set the parameters of the interaction.'
        pmi.invoke(self.worker.set,
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

    def getEpsilon(self) :
        return self.worker.getEpsilon()
    def getSigma(self) :
        return self.worker.getSigma()
    def getCutoff(self) :
        return self.worker.getCutoff()

    # define properties
    epsilon = property(getEpsilon, setEpsilon)
    sigma = property(getSigma, setSigma)
    cutoff = property(getCutoff, setCutoff)

    def computeEnergy(self, r) :
        'Compute and return the energy at the radius r.'
        return self.worker.computeEnergy(r)

    def computeForce(self, r) :
        'Compute and return the force at the radius r.'
        f = self.worker.computeForce(r)
        # TODO: Is this the right way to do it???
        f.__class__=espresso.esutil.Real3D
        return f

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
    
