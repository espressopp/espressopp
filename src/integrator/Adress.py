"""
************************************
**AdResS** - Object
************************************

The AdResS object is an extension to the integrator. It makes sure that the
integrator also processes the atomistic particles and not only the CG particles.
Hence, this object is of course only used when performing AdResS or H-AdResS
simulations.

In detail the AdResS extension makes sure:
---------------------------------------------

* that also the forces on the atomistic particles are initialized and set to
  by Adress::initForces
* that also the atomistic particles are integrated and propagated by
  Adress::integrate1 and Adress::integrate2

Example - how to turn on the AdResS integrator extension:

>>> adress      = espresso.integrator.Adress(system)
>>> integrator.addExtension(adress)

"""

from espresso.esutil import cxxinit
from espresso import pmi

from espresso.integrator.Extension import *
from _espresso import integrator_Adress

class AdressLocal(ExtensionLocal, integrator_Adress):
    'The (local) AdResS'

    def __init__(self, _system, _verletlist, _fixedtuplelist, KTI = False):
        'Local construction of a verlet list for AdResS'
        if pmi.workerIsActive():
            cxxinit(self, integrator_Adress, _system, _verletlist, _fixedtuplelist, KTI)

if pmi.isController:
    class Adress(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.integrator.AdressLocal' #,
            #pmiproperty = [ 'builds' ],
            #pmicall = [ 'totalSize', 'exclude', 'addAdrParticles', 'rebuild' ]
            )
