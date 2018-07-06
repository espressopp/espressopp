#  Copyright (C) 2012,2013,2017
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


r"""

.. function:: espressopp.analysis.Observable.compute()

		:rtype:
"""
from espressopp import pmi
from _espressopp import analysis_Observable

class result_types:
    none, real_scalar, int_scalar, real_vector, int_vector = range(5)

class ObservableLocal(object):

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

        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "compute" ]
            )
