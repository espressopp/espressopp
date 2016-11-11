#  Copyright (C) 2012,2013,2016
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
********************************************
espressopp.analysis.PressureTensorMultiLayer
********************************************

This class computes the pressure tensor of the system in `n` layers.
Layers are perpendicular to Z direction and are equidistant(distance is Lz/n).
It can be used as standalone class in python as well as
in combination with the integrator extension ExtAnalyze.

Example of standalone Usage:

>>> pt = espressopp.analysis.PressureTensorMultiLayer(system, n, dh)
>>> for i in xrange(n):
>>>     print "pressure tensor in layer %d: %s" % ( i, pt.compute())

or

>>> pt = espressopp.analysis.PressureTensorMultiLayer(system, n, dh)
>>> for k in xrange(100):
>>>     integrator.run(100)
>>>     pt.performMeasurement()
>>> for i in xrange(n):
>>>     print "average pressure tensor in layer %d: %s" % ( i, pt.compute())

Example of usage in integrator with ExtAnalyze:

>>> pt           = espressopp.analysis.PressureTensorMultiLayer(system, n, dh)
>>> extension_pt = espressopp.integrator.ExtAnalyze(pt , interval=100)
>>> integrator.addExtension(extension_pt)
>>> integrator.run(10000)
>>> pt_ave = pt.getAverageValue()
>>> for i in xrange(n):
>>>   print "average Pressure Tensor = ", pt_ave[i][:6]
>>>   print "          std deviation = ", pt_ave[i][6:]
>>> print "number of measurements  = ", pt.getNumberOfMeasurements()

The following methods are supported:

* performMeasurement()
    computes the pressure tensor and updates average and standard deviation
* reset()
    resets average and standard deviation to 0
* compute()
    computes the instant pressure tensor in `n` layers, return value: [xx, yy, zz, xy, xz, yz]
* getAverageValue()
    returns the average pressure tensor and the standard deviation,
    return value: [xx, yy, zz, xy, xz, yz, +-xx, +-yy, +-zz, +-xy, +-xz, +-yz]
* getNumberOfMeasurements()
    counts the number of measurements that have been computed (standalone or in integrator)
    does _not_ include measurements that have been done using "compute()"

.. function:: espressopp.analysis.PressureTensorMultiLayer(system, n, dh)

		:param system:
		:param n:
		:param dh:
		:type system:
		:type n:
		:type dh:
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.AnalysisBase import *
from _espressopp import analysis_PressureTensorMultiLayer

class PressureTensorMultiLayerLocal(AnalysisBaseLocal, analysis_PressureTensorMultiLayer):

    def __init__(self, system, n, dh):
	if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_PressureTensorMultiLayer, system, n, dh)

if pmi.isController:
    class PressureTensorMultiLayer(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.PressureTensorMultiLayerLocal',
            pmiproperty = [ 'n', 'dh' ]
        )
