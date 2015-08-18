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
************************************
**espressopp.analysis.RDFatomistic**
************************************


.. function:: espressopp.analysis.RDFatomistic(system, type1, type2, _span)

		:param system: 
		:param type1: 
		:param type2: 
		:param _span: 
		:type system: 
		:type type1: 
		:type type2: 
		:type _span: 

.. function:: espressopp.analysis.RDFatomistic.compute(rdfN)

		:param rdfN: 
		:type rdfN: 
		:rtype:
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_RDFatomistic

class RDFatomisticLocal(ObservableLocal, analysis_RDFatomistic):

  def __init__(self, system, type1, type2, _span):
    cxxinit(self, analysis_RDFatomistic, system, type1, type2, _span)
    
  def compute(self, rdfN):
    return self.cxxclass.compute(self, rdfN)
    
if pmi.isController :
  class RDFatomistic(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "compute" ],
      cls = 'espressopp.analysis.RDFatomisticLocal'
    )
