#  Copyright (C) 2012,2013
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
****************************************************************************
**LBOutput** - abstract base class for analysis / output in LB simulations
****************************************************************************

Abstract base class for arbitrary output from LB simulations. At the moment, the implemented realisations are:

* :class:`espressopp.analysis.LBOutputScreen` to output local density :math:`rho` and :math:`v_z` component of the velocity as a function of the coordinate :math:`x`.
* :class:`espressopp.analysis.LBOutputVzInTime` to output velocity component :math:`v_z` of a specific lattice site (the value used at the moment is :math:`{0.25*N_i, 0, 0}`) in time.
* :class:`espressopp.analysis.LBOutputVzOfX` to output simulation progress and control flux conservation when using MD to LB coupling.

.. note::

	Other types of output classes are possible. It is a subject of user requests.

..
	This happens in the plane where y and z are equal to zero (index j=k=0). The output takes place into the file vz_of_x.`step`.dat

	Computes and outputs

	* `LBOutputScreen()`
	Outputs useful information onto the screen.

	* `LBOutputVzInTime()`
	Computes and outputs a vz component of the velocity as a function of time.


.. function:: espressopp.analysis.LBOutput()

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.AnalysisBase import *
from _espressopp import analysis_LBOutput

class LBOutputLocal(AnalysisBaseLocal, analysis_LBOutput):
	#    'The (local) compute of LBOutput.'
    def writeOutput(self):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.writeOutput(self)
        
if pmi.isController :
    class LBOutput(AnalysisBase):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
#            cls =  'espressopp.analysis.LBOutputLocal',
#            pmicall = ["writeOutput"]
            )
        