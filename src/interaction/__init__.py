from espresso import pmi

from espresso.interaction.LennardJones import *

# from _espresso import interaction_FENE as _FENE
# class __FENE(FENELocal) :
#     'The FENE interaction.'
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



