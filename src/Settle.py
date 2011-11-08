from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class SettleLocal(_espresso.Settle):
    'The (local) settle.'

    def __init__(self, storage, integrator, mO=16.0, mH=1.0, distHH=1.58, distOH=1.0):
        'Local construction of a fixed touple list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.Settle, storage, integrator, mO, mH, distHH, distOH)

    def addMolecules(self, moleculelist):
        """
        Each processor takes the broadcasted list.
        """
        if pmi.workerIsActive():
            for pid in moleculelist: 
                self.cxxclass.add(self, pid)


if pmi.isController:
    class Settle(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.SettleLocal',
            pmicall = [ "addMolecules" ]
            )
