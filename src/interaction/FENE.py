from espresso import pmi
from espresso.esutil import choose

from _espresso import interaction_FENE as _FENE

class FENELocal(_FENE) :
    'The (local) FENE interaction.'
    def __init__(self, K=1.0, r0=0.0, rMax=1.0) :
        """Initialize the FENE potential.

        The parameters are identical to set()."""
        _FENE.__init__(self)
        # set the defaults
        self.set(K, r0, rMax)

    # define setter
    def set(self, K=None, r0=None, rMax=None) :
        """set((real)K, (real)r0, (real)rMax) -- Set the parameters of the interaction.
        """
        return _FENE.set(self,
            choose(K, self.K),
            choose(r0, self.r0),
            choose(rMax, self.rMax)
            )

    # define properties
    @property
    def K(self): return self.getK()
    @K.setter
    def K(self, _K): self.set(K=_K)

    @property
    def r0(self): return self.getR0()
    @r0.setter
    def r0(self, _r0): self.set(r0=_r0)

    @property
    def rMax(self): return self.getRMax()
    @K.setter
    def rMax(self, _rMax): self.set(rMax=_rMax)

# wrap FENE
if pmi.IS_CONTROLLER:
    pmi.exec_('from espresso.interaction.FENE import FENELocal')

    class FENE(object) :
        'The FENE interaction.'
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'subjectclass': 'FENELocal',
            'localcall' : [ 'computeForce', 'computeEnergy' ],
            'pmicall' : [ 'set' ],
            'pmiproperty' : ['K', 'r0', 'rMax' ]
            }


