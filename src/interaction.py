from espresso import esutil
from esutil import choose
from espresso import pmi

# wrap LennardJones
pmi.exec_('from _espresso import interaction_LennardJones as LennardJones')
class LennardJones (object):
    'The Lennard-Jones interaction.'

    def __init__(self, epsilon=1.0, sigma=1.0, cutoff=2.0) :
        """Initialize the (parallel) Lennard Jones object. 

        The parameters are identical to set."""
        object.__init__(self)
        # create the pmi object
        self.worker = pmi.create('LennardJones')
        # set the defaults
        self.set(epsilon, sigma, cutoff)

    # define setter
    def set(self, epsilon=None, sigma=None, cutoff=None) :
        """set(integer, integer, integer) -- Set the "parameters" of the interaction.
        """
        pmi.invoke('LennardJones.set',
                   self.worker,
                   choose(epsilon, self.epsilon),
                   choose(sigma, self.sigma),
                   choose(cutoff, self.cutoff)
                   )

    # define single property setters
    # avoid using these if possible
    def __setEpsilon(self, _epsilon) :
        self.set(epsilon=_epsilon)
    def __setSigma(self, _sigma) : 
        self.set(sigma=_sigma)
    def __setCutoff(self, _cutoff) : 
        self.set(cutoff=_cutoff)

    def __getEpsilon(self) :
        return self.worker.getEpsilon()
    def __getSigma(self) :
        return self.worker.getSigma()
    def __getCutoff(self) :
        return self.worker.getCutoff()

    # define properties
    epsilon = property(__getEpsilon, __setEpsilon)
    sigma = property(__getSigma, __setSigma)
    cutoff = property(__getCutoff, __setCutoff)

    def computeEnergy(self, r) :
        'Compute and return the energy at the radius r.'
        return self.worker.computeEnergy(r)

    def computeForce(self, r) :
        'Compute and return the force at the radius r.'
        f = self.worker.computeForce(r)
        return f

# wrap FENE
pmi.exec_('from _espresso import interaction_FENE as FENE')
class FENE(object) :
    'The FENE interaction.'

    def __init__(self, K=1.0, r0=0.0, rMax=1.0) :
        """Initialize the (parallel) Lennard Jones object. 

        The parameters are identical to set."""
        object.__init__(self)
        self.worker = pmi.create('FENE')
        self.set(K, r0, rMax)

    # define setter
    def set(self, K=None, r0=None, rMax=None) :
        'Set the parameters of the interaction.'
        pmi.invoke(
            'FENE.set',
            self.worker,
            choose(K, self.K),
            choose(r0, self.r0),
            choose(rMax, self.rMax)
            )

    # define single property setters
    # avoid using these if possible
    def __setK(self, _K) :
        self.set(K=_K)
    def __setR0(self, _r0) : 
        self.set(r0=_r0)
    def __setRMax(self, _rMax) : 
        self.set(rMax=_rMax)

    def __getK(self) :
        return self.worker.getK()
    def __getR0(self) :
        return self.worker.getR0()
    def __getRMax(self) :
        return self.worker.getRMax()

    # define properties
    K = property(__getK, __setK)
    r0 = property(__getR0, __setR0)
    rMax = property(__getRMax, __setRMax)

    def computeEnergy(self, r) :
        'Compute and return the energy at the radius r.'
        return self.worker.computeEnergy(r)
    
    def computeForce(self, r) :
        'Compute and return the force at the radius r.'
        f = self.worker.computeForce(r)
        # TODO: Is this the right way to do it???
        f.__class__=espresso.esutil.Real3D
        return f



