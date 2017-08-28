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

This abstract base class provides the interface and some basic
functionality for classes that do analysis or observable measurements

It provides the following methods:

.. function:: espressopp.analysis.AnalysisBase.compute()

		Computes the instant value of the observable.

		:rtype: a python list or a scalar

.. function:: espressopp.analysis.AnalysisBase.getAverageValue()

		Returns the average value for the observable and the standard deviation.

		:rtype: a python list

.. function:: espressopp.analysis.AnalysisBase.getNumberOfMeasurements()

		counts the number of measurements that have been performed (standalone or in integrator)
		does _not_ include measurements that have been done using "compute()"

		:rtype:

.. function:: espressopp.analysis.AnalysisBase.performMeasurement()

		Computes the observable and updates average and standard deviation

		:rtype:

.. function:: espressopp.analysis.AnalysisBase.reset()

		Resets average and standard deviation

		:rtype:
                    
"""

from espressopp import pmi

from espressopp.ParticleAccess import *
from _espressopp import analysis_AnalysisBase

class AnalysisBaseLocal(ParticleAccessLocal, analysis_AnalysisBase):

    def performMeasurement(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.performMeasurement(self)

    def reset(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.reset(self)

    def compute(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            res = self.cxxclass.compute(self)
            if len(res) > 1:
                return res
            else:
                return res[0]

    def getAverageValue(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getAverageValue(self)

    def getNumberOfMeasurements(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getNumberOfMeasurements(self)

if pmi.isController :
    class AnalysisBase(ParticleAccess):

        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmicall = [ "performMeasurement", "reset", "compute", "getAverageValue", "getNumberOfMeasurements" ]
            )
