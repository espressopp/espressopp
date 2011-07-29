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
