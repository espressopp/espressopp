"""
********************************
**espresso.analysis.Observable**
********************************

"""
from espresso import pmi
from _espresso import analysis_Observable

class result_types:
    none, real_scalar, int_scalar, real_vector, int_vector = range(5)

class ObservableLocal(object):
    """Abstract local base class for observables."""
    def compute(self):
        res_type = self.cxxclass.getResultType(self)
        if res_type == result_types.none:
            return
        elif res_type == result_types.real_scalar:
            return self.cxxclass.compute_real(self)
        elif res_type == result_types.int_scalar:
            return self.cxxclass.compute_int(self)
        elif res_type == result_types.real_vector:
            return self.cxxclass.compute_real_vector_python(self)
        elif res_type == result_types.int_vector:
            return self.cxxclass.compute_int_vector_python(self)
        else:
            return self.cxxclass.compute(self)

if pmi.isController :
    class Observable(object):
        """Abstract base class for observable."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "compute" ]
            )
