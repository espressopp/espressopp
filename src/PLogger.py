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
************************************
**PLogger** - Object
************************************

This module defines the parallel logger PLogger.
It can be used to switch on logging on all CPUs.

   espresso.PLogger.set(LoggerName, LoggerLevel)
   
LoggerName  : the name of the logger, if LoggerName='' than all loggers are set
LoggerLevel : possible values are 'FATAL', 'ERROR', 'WARN', 'INFO', 'TRACE', 'DEBUG'
              'DEBUG' produces most output
              'FATAL' produces least output

Example:

>>> espresso.PLogger.set('LennardJonesGeneric', 'INFO')
>>> pot = espresso.interaction.LennardJonesGeneric(1.0, 1.0, 12, 6, 1.12246)
>>> print pot.computeEnergy(1.0)
>>> espresso.PLogger.set('LennardJonesGeneric', 'ERROR')

"""

import espresso

def set(thelogger,level):
  hw = espresso.pmi.create('logging.getLogger',thelogger)
  espresso.pmi.call(hw,'setLevel',level)

