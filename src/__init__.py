# Set up logging
import logging
import os

logConfigFile="espresso_log.conf"
if os.path.exists(logConfigFile) :
    import logging.config
    logging.config.fileConfig(logConfigFile)
    log = logging.getRootLogger()
    log.info('Reading log config file %s', logConfigFile)
else :
    logging.basicConfig()
    log = logging.getLogger('root')
    log.info('Did not find log config file %s, using basic configuration.', logConfigFile)
                              
# Define the basic classes
from _espresso import Real3D
import esutil

# This injects additional methods into the Real3D class and pulls it
# into this module 
class __Real3D(Real3D) :
    """Basic 3D floating point vector as used by ESPResSo++.
    """
    __metaclass__ = esutil.ExtendBaseClass

    originit = Real3D.__init__
    def __init__(self, x=0.0, y=0.0, z=0.0) :
        return self.originit(x, y, z)

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

