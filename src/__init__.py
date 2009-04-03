# Set up logging
import logging
import os

logConfigFile="espresso_log.conf"
if os.path.exists(logConfigFile) :
    import logging.config
    logging.config.fileConfig(logConfigFile)
    log = logging.getLogger('root')
    log.info('Reading log config file %s', logConfigFile)
else :
    logging.basicConfig()
    log = logging.getLogger('root')
    log.info('Did not find log config file %s, using basic configuration.', logConfigFile)
                              
import _espresso

# This initialization routine will change existing and future loggers
# to make a connection with their Python logger and change their class

__save__logging__Logger__setLevel = logging.Logger.setLevel

def __mySetLevel__(self, level):

    __save__logging__Logger__setLevel(self, level)
    _espresso.setLogger(self)

logging.Logger.setLevel = __mySetLevel__
logging.TRACE = (logging.NOTSET + logging.DEBUG) / 2
logging.addLevelName('TRACE', logging.TRACE)

_espresso.setLogger()


# Define the basic classes
from _espresso import Real3D
import esutil

# This injects additional methods into the Real3D class and pulls it
# into this module 
class __Real3D(Real3D) :
    """Basic 3D floating point vector as used by ESPResSo++.
    """
    __metaclass__ = esutil.ExtendBaseClass

    __originit = Real3D.__init__
    def __init__(self, x=0.0, y=0.0, z=0.0) :
        return self.__originit(x, y, z)

    # create setters and getters
    @esutil.propget
    def x(self) : return self[0]

    @esutil.propset
    def x(self, v) : self[0] = v

    @esutil.propget
    def y(self) : return self[1]

    @esutil.propset
    def y(self, v) : self[1] = v

    @esutil.propget
    def z(self) : return self[2]

    @esutil.propset
    def z(self, v) : self[2] = v

    # string conversion
    def __str__(self) :
        return str(tuple(self))

