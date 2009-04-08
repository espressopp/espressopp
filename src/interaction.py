from espresso import esutil
from esutil import choose
from espresso import pmi

# Extend the C++ LennardJones class
from _espresso import interaction_LennardJones as _LennardJones

class LennardJonesLocal(_LennardJones) :
    'The (local) Lennard-Jones interaction.'

    def __init__(self, epsilon=1.0, sigma=1.0, cutoff=2.0) :
        """Initialize the local Lennard Jones object. 

        The parameters are identical to set()."""
        _LennardJones.__init__(self)
        self.set(epsilon, sigma, cutoff)

    def set(self, epsilon=None, sigma=None, cutoff=None) :
        """set( (float)epsilon, (float)sigma, (float)cutoff ) -> None -- Set the "parameters" of the interaction.
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

    # class LennardJones(object) :
    #     __metaclass__ = PMIProxy
    #     pmiclass = 'LennardJonesLocal'
    #     pmicall = [ 'set' ]
    #     pmiinvoke = []
    #     pmilocal = [ pmi.ALL ]

    import functools
    def pmi_create(cls) :
        def wrap_wrap(f) :
            @functools.wraps(f)
            def wrapper(self, *args, **keywds):
                self.local = pmi.create(cls, *args, **keywds)
                self.pmi_class = cls
                return f(self, *args, **keywds)
            return wrapper
        return wrap_wrap

#     def pmi_call(f) :
#         @functools.wraps(f)
#         def wrapper(self, *args, **kwds) :
#             pmi.call(f, self.local, *args, **kwds)
#         return wrapper
    
    pmi.exec_('from espresso.interaction import LennardJonesLocal')
    class LennardJones (object):
        'The Lennard-Jones interaction.'
        
        def __init__(self, epsilon=1.0, sigma=1.0, cutoff=2.0) :
            self.local = pmi.create('LennardJonesLocal', epsilon, sigma, cutoff)
            return object.__init__(self)
        
        def set(self, epsilon=None, sigma=None, cutoff=None) :
            pmi.call('LennardJonesLocal.set',
                     self.local, epsilon, sigma, cutoff)
            
        @property
        def epsilon(self): return self.local.epsilon
        @epsilon.setter
        def epsilon(self, _epsilon):
            pmi.call('LennardJonesLocal.epsilon.fset', self.local, _epsilon)

        @property
        def sigma(self): return self.local.sigma
        @sigma.setter
        def sigma(self, _sigma):
            pmi.call('LennardJonesLocal.sigma.fset', self.local, _sigma)

        @property
        def cutoff(self): return self.local.cutoff
        @cutoff.setter
        def cutoff(self, _cutoff):
            pmi.call('LennardJonesLocal.cutoff.fset', self.local, _cutoff)

        def computeForce(self, r) :
            return self.local.computeForce(r)

        def computeEnergy(self, r) :
            return self.local.computeEnergy(r)


# from _espresso import interaction_FENE as _FENE
# class __FENE(FENELocal) :
#     'The FENE interaction.'
#     __metaclass__ = esutil.ExtendBaseClass

#     __originit = FENELocal.__init__
#     def __init__(self, K=1.0, r0=0.0, rMax=1.0) :
#         """Initialize the FENE potential.

#         The parameters are identical to set()."""
#         self.__originit()
#         # set the defaults
#         self.set(K, r0, rMax)

#     # define setter
#     __origset = FENELocal.set
#     def set(self, K=None, r0=None, rMax=None) :
#         """set((real)K, (real)r0, (real)rMax) -- Set the "parameters" of the interaction.
#         """
#         self.__origset(
#             choose(epsilon, self.epsilon),
#             choose(sigma, self.sigma),
#             choose(cutoff, self.cutoff)
#             )

#     # define single property setters
#     # avoid using these if possible
#     def _setK(self, _K) :
#         self.set(K=_K)
#     def _setR0(self, _r0) : 
#         self.set(r0=_r0)
#     def _setRMax(self, _rMax) : 
#         self.set(rMax=_rMax)

#     # define properties
#     K = property(getK, _setK)
#     r0 = property(getR0, _setR0)
#     rMax = property(getRMax, setRMax)


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



