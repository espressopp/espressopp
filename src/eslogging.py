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

from _espresso import setLogger

# This initialization routine will change existing and future loggers
# to make a connection with their Python logger and change their class
__orig_setLevel = logging.Logger.setLevel

def __my_setLevel(self, level):
    __orig_setLevel(self, level)
    setLogger(self)

import math

logging.Logger.setLevel = __my_setLevel
logging.TRACE = int((logging.NOTSET + logging.DEBUG)/2.0)
logging.addLevelName('TRACE', logging.TRACE)
setLogger()
