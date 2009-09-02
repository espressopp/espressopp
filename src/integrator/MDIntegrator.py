from espresso import pmi

class MDIntegratorLocal(object):
    def integrate(self, steps):
        self.cxxclass.integrate(self, steps)

if pmi.IS_CONTROLLER:
    class MDIntegrator(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmiproperty = [ 'set', 
                            'velProperty', 'forceProperty', 
                            'timestep' ],
            pmicall = ['integrate']
            )
