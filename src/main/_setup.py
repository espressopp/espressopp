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


# now load the fundamental modules
# load mpi4py (must be loaded before _espressopp)

import mpi4py.MPI as MPI

# load the ES++-C++ module
import _espressopp

import logging
import logging.config
import os

# load PMI explicitly from espressopp
from espressopp import pmi

# define pmiimport
if pmi.isController:
    def pmiimport(module):
        pmi.exec_('import ' + module)
else:
    def pmiimport(module):
        pass


# set up logging
def _setupLogging():
    logConfigFile = "espressopp_log.conf"
    if os.path.exists(logConfigFile):
        logging.config.fileConfig(logConfigFile)
        log = logging.getLogger('root')
        log.info('Reading log config file %s', logConfigFile)
    else:
        logging.basicConfig(
            format="%(process)d %(asctime)s %(name)s (%(filename)s::%(lineno)s,%(funcName)s) %(levelname)s: %(message)s")
        log = logging.getLogger('root')
        log.info('Did not find log config file %s, using basic configuration.', logConfigFile)

        # This initialization routine will change existing and future loggers
        # to make a connection with their Python logger and change their class

        __orig_setLevel = logging.Logger.setLevel

        def my_setLevel(self, level):
            __orig_setLevel(self, level)
            _espressopp.setLogger(self)

        logging.Logger.setLevel = my_setLevel
        logging.TRACE = int((logging.NOTSET + logging.DEBUG) / 2.0)
        logging.addLevelName('TRACE', logging.TRACE)
        _espressopp.setLogger()


# execute the function
_setupLogging()
