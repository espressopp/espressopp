from espressopp import pmi

class SampleBaseLocal:
    def exxx(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            print('SampleBaseLocal.exxx')
            return 'abc'
    def fxxx(self):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            print('SampleBaseLocal.fxxx')
            return 'abc'

if pmi.isController :
    class SampleBase(metaclass=pmi.Proxy):
        pmiproxydefs = dict(pmicall=['exxx'])

        def fxxx(self):
            pmi.call(self.pmiobject, 'fxxx')
