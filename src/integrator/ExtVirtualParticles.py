#! /usr/bin/python
"""
Example - how to turn on the integrator extension:

>>> adress      = espressopp.integrator.ExtVirtualParticles(system)
>>> integrator.addExtension(adress)

"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_ExtVirtualParticles

class ExtVirtualParticlesLocal(ExtensionLocal, integrator_ExtVirtualParticles):
    'The (local) AdResS'

    def __init__(self, system, cl):
        'construction of a verlet list of virtual particles'
        if pmi.workerIsActive():
            cxxinit(self, integrator_ExtVirtualParticles, system, cl)

    def getCellList(self):
        'get number of interactions of the system'
        if pmi.workerIsActive():
            return self.cxxclass.getCellList(self)

    def addVirtualParticleTypes(self, pids):
        """
        Each processor takes the broadcasted atomistic particles
        and adds it to its list.
        """
        if pmi.workerIsActive():
            for pid in pids:
                self.cxxclass.addVirtualParticleType(self, pid)

if pmi.isController:
    class ExtVirtualParticles(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.integrator.ExtVirtualParticlesLocal',
            #pmiproperty = [ 'builds' ],
            pmicall = [ 'addVirtualParticleTypes', 'setFixedTupleList' ]
            )
