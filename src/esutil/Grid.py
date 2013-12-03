"""
************************
**espresso.esutil.Grid**
************************

"""
from espresso import pmi
from _espresso import esutil_Grid

class GridLocal(esutil_Grid):
  pass

if pmi.isController:
    class Grid(object):
        __metaclass__ = pmi.Proxy
        'Grid class'
        pmiproxydefs = dict(
            cls = 'espresso.esutil.GridLocal',
            localcall = [ 'mapIndexToPosition' ]
            #localcall = [ '__call__', 'normal', 'gamma', 'uniformOnSphere' ],
            #pmicall = [ 'seed' ]
        )
    
