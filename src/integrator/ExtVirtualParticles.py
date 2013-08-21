#! /usr/bin/python
"""
Example - how to turn on the integrator extension:

>>> adress      = espresso.integrator.ExtVirtualParticles(system)
>>> integrator.addExtension(adress)

"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *
from _espresso import integrator_ExtVirtualParticles

class ExtVirtualParticlesLocal(ExtensionLocal, integrator_ExtVirtualParticles):
    'The (local) AdResS'

    def __init__(self, system):
        'construction of a verlet list of virtual particles'
        if pmi.workerIsActive():
            cxxinit(self, integrator_ExtVirtualParticles, system)

if pmi.isController:
    class ExtVirtualParticles(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.integrator.ExtVirtualParticlesLocal' #,
            #pmiproperty = [ 'builds' ],
            #pmicall = [ 'totalSize', 'exclude', 'addAdrParticles', 'rebuild' ]
            )
