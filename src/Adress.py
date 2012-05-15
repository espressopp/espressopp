from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class AdressLocal(_espresso.Adress):
    'The (local) AdResS'

    def __init__(self, system, integrator):
        'Local construction of a verlet list for AdResS'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.Adress, system, integrator)

if pmi.isController:
    class Adress(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.AdressLocal' #,
            #pmiproperty = [ 'builds' ],
            #pmicall = [ 'totalSize', 'exclude', 'addAdrParticles', 'rebuild' ]
            )
