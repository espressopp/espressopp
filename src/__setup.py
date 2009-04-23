# set the path
from espresso import __setupPath

# load boostmpi (must be loaded before _espresso)
import boostmpi

# set up logging
import logging
import os
import math
import _espresso

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


# This initialization routine will change existing and future loggers
# to make a connection with their Python logger and change their class
__orig_setLevel = logging.Logger.setLevel

def __my_setLevel(self, level):
    __orig_setLevel(self, level)
    _espresso.setLogger(self)

logging.Logger.setLevel = __my_setLevel
logging.TRACE = int((logging.NOTSET + logging.DEBUG)/2.0)
logging.addLevelName('TRACE', logging.TRACE)
_espresso.setLogger()
