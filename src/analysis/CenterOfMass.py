from espresso.esutil import cxxinit
from espresso import pmi

from espresso.analysis.Observable import *
from _espresso import analysis_CenterOfMass

class CenterOfMassLocal(ObservableLocal, analysis_CenterOfMass):
    'The (local) compute of center-of-mass.'
    def __init__(self, system):
        cxxinit(self, analysis_CenterOfMass, system)

if pmi.isController :
    class CenterOfMass(Observable):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.analysis.CenterOfMassLocal'
            )
