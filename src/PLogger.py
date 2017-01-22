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
******************
espressopp.PLogger
******************

This module defines the parallel logger PLogger.
It can be used to switch on logging on all CPUs.

   espressopp.PLogger.set(LoggerName, LoggerLevel)
   
LoggerName  : the name of the logger, if LoggerName='' than all loggers are set
LoggerLevel : possible values are 'FATAL', 'ERROR', 'WARN', 'INFO', 'TRACE', 'DEBUG'
              'DEBUG' produces most output
              'FATAL' produces least output

Example:

>>> espressopp.PLogger.set('LennardJonesGeneric', 'INFO')
>>> pot = espressopp.interaction.LennardJonesGeneric(1.0, 1.0, 12, 6, 1.12246)
>>> print pot.computeEnergy(1.0)
>>> espressopp.PLogger.set('LennardJonesGeneric', 'ERROR')


.. function:: espressopp.set(thelogger, level)

		:param thelogger: 
		:param level: 
		:type thelogger: 
		:type level: 
"""

import espressopp

def set(thelogger,level):
  hw = espressopp.pmi.create('logging.getLogger',thelogger)
  espressopp.pmi.call(hw,'setLevel',level)
