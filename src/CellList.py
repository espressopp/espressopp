from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit
from math import sqrt

class CellListLocal(_espressopp.CellList):
    def __init__(self):
        'Local construction of a fixed pair list'
        if pmi.workerIsActive():
            cxxinit(self, _espressopp.CellList)

if pmi.isController:
    class CellList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.CellListLocal',
            #localcall = [ 'add' ],
            pmicall = [],
            pmiinvoke = []
        )
        
