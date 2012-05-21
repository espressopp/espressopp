from espresso import pmi
from _espresso import integrator_MDIntegrator

class MDIntegratorLocal(object):
    """Abstract local base class for molecular dynamics integrator."""
    def run(self, niter):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.run(self, niter)

    def addExtension(self, extension):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            
            # set integrator and connect to it
            extension.cxxclass.setIntegrator(extension, self)
            extension.cxxclass.connect(extension)
            
            return self.cxxclass.addExtension(self, extension)

    

if pmi.isController :
    class MDIntegrator(object):
        """Abstract base class for molecular dynamics integrator."""
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            pmiproperty = [ 'dt', 'step' ],
            pmicall = [ 'run', 'addExtension' ]
            )
