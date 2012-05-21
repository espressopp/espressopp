from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *
from _espresso import integrator_Adress

class AdressLocal(ExtensionLocal, integrator_Adress):
    'The (local) AdResS'

    def __init__(self, system):
        'Local construction of a verlet list for AdResS'
        if pmi.workerIsActive():
            cxxinit(self, integrator_Adress, system)

if pmi.isController:
    class Adress(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.integrator.AdressLocal' #,
            #pmiproperty = [ 'builds' ],
            #pmicall = [ 'totalSize', 'exclude', 'addAdrParticles', 'rebuild' ]
            )
