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

# load PMI explicitly from espressopp
from espressopp import pmi

# define pmiimport
if pmi.isController :
    def pmiimport(module):
        pmi.exec_('import ' + module)
else:
    def pmiimport(module):
        pass

# set up logging
def _setupLogging():
    import logging, os, math

    logConfigFile="espressopp_log.conf"
    if os.path.exists(logConfigFile) :
        import logging.config
        logging.config.fileConfig(logConfigFile)
        log = logging.getLogger('root')
        log.info('Reading log config file %s', logConfigFile)
    else :
        logging.basicConfig(
           format = "%(process)d %(asctime)s %(name)s (%(filename)s::%(lineno)s,%(funcName)s) %(levelname)s: %(message)s")
        log = logging.getLogger('root')
        log.info('Did not find log config file %s, using basic configuration.', logConfigFile)

        # This initialization routine will change existing and future loggers
        # to make a connection with their Python logger and change their class

        def __my_setLevel(self, level):
            __orig_setLevel(self, level)
            _espressopp.setLogger(self)

        __orig_setLevel = logging.Logger.setLevel
        logging.Logger.setLevel = __my_setLevel
        logging.TRACE = int((logging.NOTSET + logging.DEBUG)/2.0)
        logging.addLevelName('TRACE', logging.TRACE)
        _espressopp.setLogger()

# execute the function

_setupLogging()

def _setupProperty():
    import __builtin__

    # Make the property setter decorator syntax of python 2.6+ available
    # to earlier versions
    try :
        __setter = __builtin__.property.setter
    except AttributeError :
        import __builtin__, sys
        # save the property builtin
        _property = __builtin__.property
        # now define our property
        # stolen from http://bruynooghe.blogspot.com/2008/04/xsetter-syntax-in-python-25.html 
        class property(_property):
            def __init__(self, fget, *args, **kwargs):
                self.__doc__ = fget.__doc__
                super(property, self).__init__(fget, *args, **kwargs)

            def setter(self, fset):
                cls_ns = sys._getframe(1).f_locals
                for k, v in cls_ns.iteritems():
                    if v == self:
                        propname = k
                        break
                cls_ns[propname] = _property(self.fget, fset,
                                            self.fdel, self.__doc__)
                return cls_ns[propname]

            def deleter(self, fdel):
                cls_ns = sys._getframe(1).f_locals
                for k, v in cls_ns.iteritems():
                    if v == self:
                        propname = k
                        break
                cls_ns[propname] = _property(self.fget, self.fset,
                                            fdel, self.__doc__)
                return cls_ns[propname]

        # Now override the property builtin
        __builtin__.property = property

_setupProperty()
