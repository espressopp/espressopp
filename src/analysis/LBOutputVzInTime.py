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
****************************************
**espressopp.analysis.LBOutputVzInTime**
****************************************
Child class derived from the abstract class 'LBOutput()'. For more information and examples of usage, please, refer to this class.
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.LBOutput import *
from _espressopp import analysis_LBOutput_VzInTime

class LBOutputVzInTimeLocal(LBOutputLocal, analysis_LBOutput_VzInTime):
	#    'The (local) compute of LBOutputVzInTime.'
    def __init__(self, system, latticeboltzmann):
        if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, analysis_LBOutput_VzInTime, system, latticeboltzmann)

if pmi.isController :
    class LBOutputVzInTime(LBOutput):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.analysis.LBOutputVzInTimeLocal',
            pmicall = ["writeOutput"]
            )