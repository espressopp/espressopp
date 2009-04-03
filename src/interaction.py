from espresso import esutil
from esutil import choose
from espresso import pmi
from _espresso import interaction_LennardJones as LennardJonesLocal

class __LennardJones(LennardJonesLocal) :
    'The Lennard-Jones interaction.'
    __metaclass__ = esutil.ExtendBaseClass

    __originit = LennardJonesLocal.__init__
    def __init__(self, epsilon=1.0, sigma=1.0, cutoff=2.0) :
        """Initialize the (parallel) Lennard Jones object. 

        The parameters are identical to set."""
        self.__originit()
        # set the defaults
        self.set(epsilon, sigma, cutoff)

    # define setter
    __origset = LennardJonesLocal.set
    def set(self, epsilon=None, sigma=None, cutoff=None) :
        """set(integer, integer, integer) -- Set the "parameters" of the interaction.
        """
        self.__origset(
            choose(epsilon, self.epsilon),
            choose(sigma, self.sigma),
            choose(cutoff, self.cutoff)
            )

    # define single property setters
    # avoid using these if possible
    def _setEpsilon(self, _epsilon) :
        self.set(epsilon=_epsilon)
    def _setSigma(self, _sigma) : 
        self.set(sigma=_sigma)
    def _setCutoff(self, _cutoff) : 
        self.set(cutoff=_cutoff)

    def _getEpsilon(self) :
        return self.getEpsilon()
    def _getSigma(self) :
        return self.getSigma()
    def _getCutoff(self) :
        return self.getCutoff()

    # define properties
    epsilon = property(_getEpsilon, _setEpsilon)
    sigma = property(_getSigma, _setSigma)
    cutoff = property(_getCutoff, _setCutoff)


if pmi.IS_CONTROLLER :
    # wrap LennardJones

    # class LennardJones(object) :
    #     __metaclass__ = PMIProxy
    #     pmiclass = 'LennardJonesLocal'
    #     pmicall = [ 'set' ]
    #     pmiinvoke = []
    #     pmilocal = [ pmi.ALL ]
    
    pmi.exec_('from espresso.interaction import LennardJonesLocal')
    class LennardJones (object):
        'The Lennard-Jones interaction.'
        
        def __init__(self, epsilon=1.0, sigma=1.0, cutoff=2.0) :
            """Initialize the (parallel) Lennard Jones object. 
            
            The parameters are identical to set."""
            # create the pmi object
            self.local = pmi.create('LennardJonesLocal',
                                    epsilon, sigma, cutoff)
            return object.__init__(self)
        
        # define setter
        def set(self, epsilon=None, sigma=None, cutoff=None) :
            """set(integer, integer, integer) -- Set the "parameters" of the interaction.
            """
            pmi.call('LennardJonesLocal.set',
                     self.local, epsilon, sigma, cutoff)
            
        # define single property setters
        # avoid using these if possible
        def _setEpsilon(self, _epsilon) :
            pmi.call('LennardJonesLocal._setEpsilon',
                     self.local, _epsilon)
        def _setSigma(self, _sigma) : 
            pmi.call('LennardJonesLocal._setSigma',
                     self.local, _sigma)
        def _setCutoff(self, _cutoff) : 
            pmi.call('LennardJonesLocal._setCutoff',
                     self.local, _cutoff)
                        
        def _getEpsilon(self) :
            return self.local.getEpsilon()
        def _getSigma(self) :
            return self.local.getSigma()
        def _getCutoff(self) :
            return self.local.getCutoff()
                        
        # define properties
        epsilon = property(_getEpsilon, _setEpsilon)
        sigma = property(_getSigma, _setSigma)
        cutoff = property(_getCutoff, _setCutoff)
                        
        def computeEnergy(self, r) :
            'Compute and return the energy at the radius r.'
            return self.local.computeEnergy(r)
                        
        def computeForce(self, r) :
            'Compute and return the force at the radius r.'
            f = self.local.computeForce(r)
            return f



# # wrap FENE
# pmi.exec_('from _espresso import interaction_FENE as FENE')
# class FENE(object) :
#     'The FENE interaction.'

#     def __init__(self, K=1.0, r0=0.0, rMax=1.0) :
#         """Initialize the (parallel) Lennard Jones object. 

#         The parameters are identical to set."""
#         object.__init__(self)
#         self.local = pmi.create('FENE')
#         self.set(K, r0, rMax)

#     # define setter
#     def set(self, K=None, r0=None, rMax=None) :
#         'Set the parameters of the interaction.'
#         pmi.invoke(
#             'FENE.set',
#             self.local,
#             choose(K, self.K),
#             choose(r0, self.r0),
#             choose(rMax, self.rMax)
#             )

#     # define single property setters
#     # avoid using these if possible
#     def __setK(self, _K) :
#         self.set(K=_K)
#     def __setR0(self, _r0) : 
#         self.set(r0=_r0)
#     def __setRMax(self, _rMax) : 
#         self.set(rMax=_rMax)

#     def __getK(self) :
#         return self.local.getK()
#     def __getR0(self) :
#         return self.local.getR0()
#     def __getRMax(self) :
#         return self.local.getRMax()

#     # define properties
#     K = property(__getK, __setK)
#     r0 = property(__getR0, __setR0)
#     rMax = property(__getRMax, __setRMax)

#     def computeEnergy(self, r) :
#         'Compute and return the energy at the radius r.'
#         return self.local.computeEnergy(r)
    
#     def computeForce(self, r) :
#         'Compute and return the force at the radius r.'
#         f = self.local.computeForce(r)
#         # TODO: Is this the right way to do it???
#         f.__class__=espresso.esutil.Real3D
#         return f



