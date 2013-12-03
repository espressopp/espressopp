"""
***********************
**espresso.Exceptions**
***********************

"""
import sys, traceback

class Error(Exception):
    """Raised to show unrecoverable espresso errors.
    """
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
    """ Raised to indicate, that a certain Particle does not exist on a CPU
    """
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
    """ Raised to indicate, that a certain Particle property does not exists 
    """
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
    """ Raised to indicate, that a FixedPairList object is missing
    """
    def __init__(self, msg):
        try:
            raise Exception
        except:
            self.msg = msg
    def __str__(self) :
        return self.msg
    def __repr__(self) :
        return str(self)
