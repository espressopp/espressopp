# set the path
from espresso import __setupPath

# load boostmpi (must be loaded before _espresso)
from espresso import boostmpi

# set up logging
import logging, os, math
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

# Make the property setter decorator syntax of python 2.6+ available
# to earlier versions
try :
    __setter = property.setter
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

