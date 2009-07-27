from espresso import pmi

class MDIntegratorLocal(object):
    pass

if pmi.IS_CONTROLLER:
    pmi.exec_('import espresso.integrator')
    class MDIntegrator(object):
        pass
#         __metaclass__ = pmi.Proxy
#         pmiproxydefs = {
#             'pmicall': ['run', 'step'],
#             'pmiproperty': [ 'timestep' ]
#             }    
