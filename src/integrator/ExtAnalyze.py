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


"""
*************************************
**ExtAnalyze** - Integrator Extension
*************************************

This class can be used to execute nearly all analysis objects
within the main integration loop which allows to automatically
accumulate time averages (with standard deviation error bars). 
  
Example Usage:
-----------------

>>> pt           = espressopp.analysis.PressureTensor(system)
>>> extension_pt = espressopp.integrator.ExtAnalyze(pt , interval=100)
>>> integrator.addExtension(extension_pt)
>>> integrator.run(10000)
>>>
>>> pt_ave = pt.getAverageValue()
>>> print "average Pressure Tensor = ", pt_ave[:6]
>>> print "          std deviation = ", pt_ave[6:]
>>> print "number of measurements  = ", pt.getNumberOfMeasurements()
"""

from espressopp.esutil import cxxinit
from espressopp import pmi
from espressopp.integrator.Extension import *
from _espressopp import integrator_ExtAnalyze 

class ExtAnalyzeLocal(ExtensionLocal, integrator_ExtAnalyze):
    'The (local) extension analyze.'
    def __init__(self, action_obj, interval=1):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
           cxxinit(self, integrator_ExtAnalyze, action_obj, interval)

if pmi.isController :
    class ExtAnalyze(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.ExtAnalyzeLocal',
        )
