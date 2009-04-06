import unittest
import logging
import logging.handlers
import os
import espresso
from espresso.interaction import LennardJones

filename = "log.out"

class Test0Logging(unittest.TestCase) :
   def setUp(self) :
      global filename
      if os.path.exists(filename):
         os.remove(filename)
      
      # create logger
      self.log = logging.getLogger("_espresso.interaction.LennardJones")
      self.log.setLevel(logging.TRACE)
      # create handler 
      handler = logging.FileHandler(filename)
      # create formatter
      formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
      # add formatter to handler
      handler.setFormatter(formatter)
      # add handler to logger
      self.log.addHandler(handler)
      
   def test0Create(self) :
      
      lj = LennardJones(1.0, 2.0, 3.0)
      self.assertEqual(lj.epsilon, 1.0)
      self.assertEqual(lj.sigma, 2.0)
      self.assertEqual(lj.cutoff, 3.0)
      
      # now read the file log.out and find "DEBUG" and "interaction"
      
      f = open(filename, "r")
      s = f.read()
      self.assert_(s.find("_espresso.interaction.LennardJones") > 0)
      self.assert_(s.find("DEBUG") > 0)
      self.assert_(s.find("TRACE") < 0)
      f.close()
   

if __name__ == "__main__":
    unittest.main()
