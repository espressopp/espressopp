from espresso.esutil import choose
from espresso import pmi

from _espresso import interaction_LennardJones as _LennardJones

class LennardJonesLocal(_LennardJones) :
    'The (local) Lennard-Jones interaction.'

    def __init__(self, epsilon=1.0, sigma=1.0, cutoff=2.0) :
        """Initialize the local Lennard Jones object. 

        The parameters are identical to set()."""
        _LennardJones.__init__(self)
        self.set(epsilon, sigma, cutoff)

    def set(self, epsilon=None, sigma=None, cutoff=None) :
        """set( (float)epsilon, (float)sigma, (float)cutoff ) -> None -- Set the \"parameters\" of the interaction.
        """
        return _LennardJones.set(self,
                                 choose(epsilon, self.epsilon),
                                 choose(sigma, self.sigma),
                                 choose(cutoff, self.cutoff)
                                 )

    # define properties
    @property
    def epsilon(self) : return self.getEpsilon()
    @epsilon.setter
    def epsilon(self, _epsilon) : self.set(epsilon=_epsilon)

    @property
    def sigma(self) : return self.getSigma()
    @sigma.setter
    def sigma(self, _sigma) : self.set(sigma=_sigma)

    @property
    def cutoff(self) : return self.getCutoff()
    @cutoff.setter
    def cutoff(self, _cutoff) : self.set(cutoff=_cutoff)

if pmi.IS_CONTROLLER :
    # wrap LennardJones
    
    pmi.exec_('from espresso.interaction.LennardJones import LennardJonesLocal')
    class LennardJones (object):
        'The Lennard-Jones interaction.'

        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'subjectclass': 'LennardJonesLocal',
            'localcall' : [ 'computeForce', 'computeEnergy' ],
            'pmicall' : [ 'set' ],
            'pmiproperty' : [ 'epsilon', 'sigma', 'cutoff' ]
            }
