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


import unittest
import logging
import logging.handlers
import os
import espressopp
from espressopp.interaction import LennardJones

filename = "log.out"

class Test0Logging(unittest.TestCase) :
   def test0Create(self) :
      
      lj = LennardJones(1.0, 2.0, 3.0)
      self.assertEqual(lj.epsilon, 1.0)
      self.assertEqual(lj.sigma, 2.0)
      self.assertEqual(lj.cutoff, 3.0)
      
      # now read the file log.out and find "DEBUG" and "potential"
      
      f = open(filename, "r")
      s = f.read()
      self.assert_(s.find("_espressopp.interaction.LennardJones") > 0)
      self.assert_(s.find("DEBUG") > 0)
      self.assert_(s.find("TRACE") < 0)
      f.close()
   

if __name__ == "__main__":
   if os.path.exists(filename):
      os.remove(filename)
      
   # create logger
   log = logging.getLogger("_espressopp.interaction.LennardJones")
   log.setLevel(logging.TRACE)
   # deactivate propagation of log messages up the hierarchy
   log.propagate=0
   # create handler 
   handler = logging.FileHandler(filename)
   handler.setLevel(logging.TRACE)
   # create formatter
   formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
   # add formatter to handler
   handler.setFormatter(formatter)
   # add handler to logger
   log.addHandler(handler)
   
   unittest.main()
