"""
*********************************
**espresso.integrator.Extension**
*********************************

"""
#from espresso.esutil import cxxinit
from espresso import pmi
from _espresso import integrator_Extension 

class ExtensionLocal(object):
    'The (local) Extension abstract base class.'
    
    #def __init__(self, integrator):
    #    if pmi.workerIsActive():    
    #        cxxinit(self, integrator)
            
    #        # set center of TD force
    #        if (center != []):
    #            self.cxxclass.setCenter(self, center[0], center[1], center[2])

    #def addForce(self, itype, filename, type):
    #        """
    #        Each processor takes the broadcasted interpolation type,
    #        filename and particle type
    #        """
    #        if pmi.workerIsActive():
    #            self.cxxclass.addForce(self, itype, filename, type)
    
    def connect(self):
      return self.cxxclass.connect(self)
    def disconnect(self):
      return self.cxxclass.disconnect(self)

if pmi.isController :
    class Extension(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            #cls =  'espresso.integrator.Extension',
            #pmiproperty = [ 'itype', 'filename'],
            #pmicall = ['addForce']
            pmicall = [ 'connect', 'disconnect' ]
        )
