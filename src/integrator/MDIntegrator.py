from espresso import pmi

from _espresso import integrator_MDIntegrator
class MDIntegratorLocal(integrator_MDIntegrator):
    pass

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.integrator')
    class MDIntegrator(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = {
            'pmicall': ['run', 'step'],
            'pmiproperty': [ 'timestep' ]
            }


    
