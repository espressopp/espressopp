from espresso import pmi
from _espresso import analysis_Observable

class ObservableLocal(object):
    """Abstract local base class for observables."""
    def compute(self):
        return self.cxxclass.compute(self)

if pmi.isController :
    class Observable(object):
        """Abstract base class for observable."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "compute" ]
            )
