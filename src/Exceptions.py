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
*********************
espressopp.Exceptions
*********************


.. function:: espressopp.Error(msg)

		:param msg: 
		:type msg: 

.. function:: espressopp.ParticleDoesNotExistHere(msg)

		:param msg: 
		:type msg: 

.. function:: espressopp.UnknownParticleProperty(msg)

		:param msg: 
		:type msg: 

.. function:: espressopp.MissingFixedPairList(msg)

		:param msg: 
		:type msg: 
"""
import sys, traceback

class Error(Exception):


    def __init__(self, msg):
        try:
            raise Exception
        except:
            file, lineno, module, line = traceback.extract_stack()[0]
            self.msg = 'ERROR while executing ' + str(file) + ' line ' + str(lineno) + ': ' + str(line) + '\n-> ' + msg
    def __str__(self) :
        return self.msg
    def __repr__(self) :
        return str(self)

class ParticleDoesNotExistHere(Exception):


    def __init__(self, msg):
        try:
            raise Exception
        except:
            self.msg = msg
    def __str__(self) :
        return self.msg
    def __repr__(self) :
        return str(self)

class UnknownParticleProperty(Exception):


    def __init__(self, msg):
        try:
            raise Exception
        except:
            self.msg = msg
    def __str__(self) :
        return self.msg
    def __repr__(self) :
        return str(self)

class MissingFixedPairList(Exception):


    def __init__(self, msg):
        try:
            raise Exception
        except:
            self.msg = msg
    def __str__(self) :
        return self.msg
    def __repr__(self) :
        return str(self)
