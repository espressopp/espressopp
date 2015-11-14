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
***********************************
**espressopp.analysis.StaticStructF**
***********************************

"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.analysis.Observable import *
from _espressopp import analysis_StaticStructF

class StaticStructFLocal(ObservableLocal, analysis_StaticStructF):
  'The (local) compute the static structure function.'
  def __init__(self, system):
    cxxinit(self, analysis_StaticStructF, system)
    
  def compute(self, nqx, nqy, nqz, bin_factor, ofile = None):
    if ofile is None:
      return self.cxxclass.compute(self, nqx, nqy, nqz, bin_factor)
    else:    
      #run compute on each CPU
      result = self.cxxclass.compute(self, nqx, nqy, nqz, bin_factor)
      #create the outfile only on CPU 0
      if pmi.isController:
        myofile = 'qsq_' + str(ofile) + '.txt'
        outfile = open (myofile, 'w')
        for i in range (len(result)):
          line = str(result[i][0]) + "\t" + str(result[i][1]) + "\n"
          outfile.write(line)
        outfile.close()
      return result  

  def computeSingleChain(self, nqx, nqy, nqz, bin_factor, chainlength, ofile = None):
    if ofile is None:
      return self.cxxclass.computeSingleChain(self, nqx, nqy, nqz, bin_factor, chainlength)
    else:
      #run computeSingleChain on each CPU
      result = self.cxxclass.computeSingleChain(self, nqx, nqy, nqz, bin_factor, chainlength)
      print result #this line is in case the outfile causes problems
      #create the outfile only on CPU 0
      if pmi.isController:
        myofile = 'qsq_singleChain' + str(ofile) + '.txt'
        outfile = open (myofile, 'w')
        for i in range (len(result)):
          line = str(result[i][0]) + "\t" + str(result[i][1]) + "\n"
          outfile.write(line)
        outfile.close()
      return result
   
if pmi.isController:
  class StaticStructF(Observable):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      pmicall = [ "compute", "computeSingleChain" ],
      cls = 'espressopp.analysis.StaticStructFLocal'
    )