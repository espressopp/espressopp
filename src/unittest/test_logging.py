import unittest

import logging
import os

filename = "log.out"

if os.path.exists(filename):
   os.remove(filename)

if len(logging.root.handlers) == 0:
   # do not change the logging level
   logging.basicConfig(format = 
     "%(asctime)s %(name)s (%(filename)s::%(lineno)s,%(funcName)s) %(levelname)s: %(message) s",
     filename = filename)

import espresso

LennardJones = espresso._espresso.interaction_LennardJones 

class Test0Logging(unittest.TestCase) :

    def test0Create(self) :

        logging.getLogger("interaction.LennardJones").setLevel(logging.TRACE)

        lj = LennardJones()
        lj.set(1.0, 2.0, 3.0)
        self.assertEqual(lj.getEpsilon(), 1.0)
        self.assertEqual(lj.getSigma(), 2.0)
        self.assertEqual(lj.getCutoff(), 3.0)

        # now read the file log.out and find "DEBUG" and "interaction"

        f = open(filename, "r")
        str = f.read()
        self.assert_(str.find("interaction") > 0)
        self.assert_(str.find("DEBUG") > 0)
        self.assert_(str.find("TRACE") < 0)
        f.close()


if __name__ == "__main__":
    unittest.main()
