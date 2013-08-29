from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit
from math import sqrt

class CellListLocal(_espresso.CellList):
    def __init__(self):
        'Local construction of a fixed pair list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.CellList)

if pmi.isController:
    class CellList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.CellListLocal',
            #localcall = [ 'add' ],
            pmicall = [],
            pmiinvoke = []
        )
        
